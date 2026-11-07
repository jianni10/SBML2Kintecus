#!/usr/bin/env python3
"""
kintcheck.py -- run every .kintpack through Kintecus.exe and compare the output with the SBML reference
(libRoadRunner) time course shipped in references/<ID>.csv.

Pure Python 3.8+ standard library; no packages to install.

    python kintcheck.py --kintecus "C:\\Kintecus\\Kintecus.exe"            (all packs under this folder)
    python kintcheck.py --kintecus ... --only BIOMD0000000010 BIOMD0000000236
    python kintcheck.py --kintecus ... --folder 1_kintecus_verified --stall 120 --cap 1800

For each pack it writes MODEL.DAT / SPECIES.DAT / PARM.DAT and any profile (O_) files exactly as the Kintecus
GUI does, runs Kintecus with the pack's switches plus -i (interactive mode: no key presses needed; KDONE.OUT /
KERR.OUT report completion or a fatal error, CONC.TXT growth is the progress signal), then scores

    match error = max over compared species and rows of |Kintecus - reference| / (that species' peak in the reference)

SUCCESS = match error under 4% (change with --tol; 0.001 = the earlier 0.1% rule). Results go to kintcheck_report.csv and kintcheck_report.txt
(written as it goes, so a stopped run keeps what it finished). The column "result_when_built" shows what the
same pack gave when it was built, so changes in a new Kintecus version stand out.
"""
import argparse, csv, glob, json, math, os, shutil, subprocess, sys, time, bisect, datetime

HERE = os.path.dirname(os.path.abspath(__file__))


# ------------------------------------------------------------------ pack -> DAT files (Kintecus GUI/VBA rules)
def write_grid(path, rows, ncols, trailing_tab):
    out = []
    for r in rows:
        cells = (list(r) + [""] * ncols)[:ncols]
        if cells[0] == "":
            cells[0] = "# "
        out.append("\t".join(cells) + ("\t" if trailing_tab else ""))
    with open(path, "w", newline="", encoding="latin-1", errors="replace") as fh:
        fh.write("\r\n".join(out) + "\r\n")


def write_deck(pack, d):
    g = pack["grids"]
    write_grid(os.path.join(d, "MODEL.DAT"), g["model"]["rows"], 5, False)
    write_grid(os.path.join(d, "SPECIES.DAT"), g["species"]["rows"], 9, True)
    write_grid(os.path.join(d, "PARM.DAT"), g["parm"]["rows"], 6, True)
    if g.get("therm", {}).get("rows"):
        write_grid(os.path.join(d, "THERM.DAT"), g["therm"]["rows"], 6, True)
    for of in pack.get("ofiles", []):
        with open(os.path.join(d, of["filename"]), "w", newline="") as fh:
            fh.write("\r\n".join("\t".join(r) for r in of["rows"]) + "\r\n")


# ------------------------------------------------------------------ reading files
def read_conc(path):
    """CONC.TXT -> (TIME list, {COLUMN_NAME_UPPER: list})."""
    with open(path, encoding="latin-1", errors="replace") as fh:
        lines = [l.rstrip("\r\n").split("\t") for l in fh if l.strip()]
    if not lines:
        return [], {}
    hdr = [h.strip() for h in lines[0]]
    cols = {h.upper(): [] for h in hdr if h}
    idx = [(i, h.upper()) for i, h in enumerate(hdr) if h]
    t = []
    for l in lines[1:]:
        try:
            vals = {name: float(l[i]) for i, name in idx if i < len(l) and l[i].strip()}
        except ValueError:
            continue
        if len(vals) < len(idx):
            continue
        for name, v in vals.items():
            cols[name].append(v)
    tname = idx[0][1]
    t = cols.pop(tname, [])
    return t, cols


def conc_progress(path):
    try:
        with open(path, "rb") as fh:
            data = fh.read()
    except OSError:
        return 0, 0.0
    lines = [l for l in data.split(b"\n") if l.strip()]
    for l in reversed(lines[1:]):
        try:
            return len(lines) - 1, float(l.split(b"\t")[0])
        except ValueError:
            continue
    return max(len(lines) - 1, 0), 0.0


def read_reference(path):
    meta, scale, names, rows = {}, [], [], []
    with open(path, newline="") as fh:
        for r in csv.reader(fh):
            if not r:
                continue
            if r[0].startswith("#"):
                key = r[0][1:]
                if key == "scale":
                    scale = [float(x) for x in r[1:]]
                else:
                    meta[key] = r[1:]
                continue
            if not names:
                names = r[1:]
                continue
            rows.append([float(x) for x in r])
    return meta, names, scale, rows


# ------------------------------------------------------------------ running Kintecus
PROGRESS = 0                 # seconds between progress lines (0 = none); set by --progress
PROGRESS_STATE = {"last": 0.0}
NO_PROJECTION = False         # set by --no-projection


def run_kintecus(exe, switches, d, runner, t_end, stall, cap, show):
    """Returns (state, message, wall_seconds). state: done | fatal | stalled | too_slow | timeout | launch_error"""
    for f in ("CONC.TXT", "KDONE.OUT", "KERR.OUT", "KRUN.TMP"):
        try:
            os.remove(os.path.join(d, f))
        except OSError:
            pass
    args = [exe] + switches.split() + ["-i"]
    kw = dict(cwd=d, stdin=subprocess.DEVNULL)
    if runner:                                     # e.g. testing under Wine: --runner "wine {cmd}"
        cmd = runner.replace("{cmd}", subprocess.list2cmdline(args) if os.name == "nt" else " ".join(args))
        kw.update(shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if os.name != "nt":
            kw["start_new_session"] = True
        target = cmd
    else:
        target = args
        if os.name == "nt":                        # own (hidden) console: Kintecus writes to the console directly
            kw["creationflags"] = subprocess.CREATE_NEW_CONSOLE
            if not show:
                si = subprocess.STARTUPINFO()
                si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                si.wShowWindow = 0                 # SW_HIDE
                kw["startupinfo"] = si
        else:
            kw.update(stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    start = time.time()
    try:
        p = subprocess.Popen(target, **kw)
    except OSError as ex:
        return "launch_error", str(ex), 0.0
    last_rows, last_change, state, msg = 0, start, None, ""
    PROGRESS_STATE["last"] = start
    conc = os.path.join(d, "CONC.TXT")
    while True:
        time.sleep(0.5)
        now = time.time()
        el = now - start
        if os.path.exists(os.path.join(d, "KERR.OUT")):
            time.sleep(1.0)
            try:
                with open(os.path.join(d, "KERR.OUT"), encoding="latin-1", errors="replace") as fh:
                    msg = " ".join(fh.read().split())
            except OSError:
                msg = "FATAL"
            state = "fatal"
            break
        if os.path.exists(os.path.join(d, "KDONE.OUT")):
            time.sleep(0.5)
            try:
                with open(os.path.join(d, "KDONE.OUT"), encoding="latin-1", errors="replace") as fh:
                    txt = " ".join(fh.read().split())
            except OSError:
                txt = ""
            state, msg = ("fatal", txt) if "FATAL" in txt.upper() else ("done", txt)
            break
        if p.poll() is not None:
            time.sleep(0.5)
            if os.path.exists(os.path.join(d, "KDONE.OUT")) or os.path.exists(os.path.join(d, "KERR.OUT")):
                continue
            state, msg = "done", f"Kintecus exited (code {p.returncode})"
            break
        rows, tsim = conc_progress(conc)
        if rows > last_rows:
            last_rows, last_change = rows, now
        if PROGRESS and now - PROGRESS_STATE["last"] >= PROGRESS:          # long jobs: show they are alive
            PROGRESS_STATE["last"] = now
            pct = f" ({100 * tsim / t_end:.1f}%)" if t_end else ""
            eta = (f", about {el * (t_end - tsim) / tsim / 60:.0f} min to go" if t_end and tsim > 0 else "")
            print(f"      ... {el / 60:6.1f} min: t = {tsim:g} of {t_end:g}{pct}{eta}", flush=True)
        if stall and now - last_change > stall and el > stall:
            state, msg = "stalled", f"no new output row for {stall:.0f} s (simulated time {tsim:g} of {t_end:g})"
        elif cap and t_end and not NO_PROJECTION and el > max(60, 0.2 * cap) and tsim > 0 and el * t_end / tsim > cap:
            state, msg = "too_slow", f"projected {el * t_end / tsim:.0f} s > cap {cap:.0f} s (at t={tsim:g} of {t_end:g})"
        elif cap and el > cap:
            state, msg = "timeout", f"over {cap:.0f} s (simulated time {tsim:g} of {t_end:g})"
        if state:
            break
    wall = time.time() - start
    if p.poll() is None:
        try:
            if runner and os.name != "nt":
                import signal
                os.killpg(p.pid, signal.SIGKILL)
            else:
                p.kill()
            p.wait(timeout=10)
        except Exception:
            pass
    return state, msg, wall


# ------------------------------------------------------------------ comparison
def compare(t_kin, cols, ref_names, scale, ref_rows, samp):
    """max over species of max |kin - ref| / scale, using rows the Kintecus run reached."""
    if not t_kin:
        return None, None, "Kintecus wrote no output rows", 0
    order = sorted(range(len(t_kin)), key=lambda i: t_kin[i])
    tk = [t_kin[i] for i in order]
    tol_t = 1e-6 * max(samp, 1e-300)
    worst, worst_sp, missing, n = 0.0, None, [], 0
    per = {}
    for j, nm in enumerate(ref_names):
        col = cols.get(nm.upper())
        if col is None:
            missing.append(nm)
            continue
        ck = [col[i] for i in order]
        e = 0.0
        for r in ref_rows:
            t, ref = r[0], r[j + 1]
            if t > tk[-1] + tol_t:
                continue
            k = bisect.bisect_left(tk, t - tol_t)
            if k < len(tk) and abs(tk[k] - t) <= tol_t:
                v = ck[k]
            else:                                  # off-grid: linear interpolation
                k = min(max(k, 1), len(tk) - 1)
                t0, t1 = tk[k - 1], tk[k]
                v = ck[k - 1] + (ck[k] - ck[k - 1]) * ((t - t0) / (t1 - t0) if t1 > t0 else 0.0)
            if not math.isfinite(v):
                e = float("inf"); break
            e = max(e, abs(v - ref) / scale[j])
            n += 1
        per[nm] = e
        if worst_sp is None or e > worst:
            worst, worst_sp = e, nm
    note = f"species missing from CONC.TXT: {', '.join(missing[:5])}" if missing else ""
    if worst_sp is None:
        return None, None, note or "no species to compare", 0
    return worst, worst_sp, note, n


# ------------------------------------------------------------------ main
def trim_big(d, limit=20_000_000, keep=1000):
    """Kintecus's screen log VIEW.TXT (and KWARN.OUT) can reach gigabytes when a warning repeats every step; they
    are not used here. Files over 20 MB in a kept run folder are cut to their first and last 1,000 lines
    (<name>_head_tail.txt) and the original is deleted."""
    for name in os.listdir(d):
        if name.upper() not in ("VIEW.TXT", "KWARN.OUT"):
            continue
        f = os.path.join(d, name)
        try:
            if os.path.getsize(f) <= limit:
                continue
            from collections import deque
            head, tail, n = [], deque(maxlen=keep), 0
            with open(f, encoding="latin-1", errors="replace") as fh:
                for line in fh:
                    n += 1
                    if n <= keep:
                        head.append(line)
                    else:
                        tail.append(line)
            base = os.path.splitext(name)[0]
            with open(os.path.join(d, base + "_head_tail.txt"), "w", encoding="latin-1") as out:
                out.writelines(head)
                out.write(f"\n... {n - len(head) - len(tail):,} lines removed by kintcheck ({os.path.getsize(f):,} bytes) ...\n\n")
                out.writelines(tail)
            os.remove(f)
        except OSError:
            pass


def main():
    ap = argparse.ArgumentParser(description="Run .kintpack files through Kintecus and compare with SBML references")
    ap.add_argument("--kintecus", default=os.environ.get("KINTECUS_EXE", os.path.join(HERE, "Kintecus.exe")),
                    help="path to Kintecus.exe")
    ap.add_argument("--packs", default=HERE, help="folder searched (recursively) for .kintpack files")
    ap.add_argument("--refs", default=os.path.join(HERE, "references"), help="folder with <ID>.csv references")
    ap.add_argument("--work", default=os.path.join(HERE, "kintcheck_work"), help="scratch folder for the runs")
    ap.add_argument("--out", default=HERE, help="where kintcheck_report.csv/.txt are written")
    ap.add_argument("--only", nargs="*", help="only these model IDs")
    ap.add_argument("--folder", help="only packs whose folder name contains this text")
    ap.add_argument("--tol", type=float, default=4e-2, help="success threshold (fraction of peak); default 0.04 (4%%)")
    ap.add_argument("--stall", type=float, default=60,
                    help="give up after this many seconds without a new row (0 = never)")
    ap.add_argument("--cap", type=float, default=900, help="wall-clock limit per pack in seconds (0 = no limit)")
    ap.add_argument("--jobs", help="CSV with columns id,switches: run only these packs, each with its own switches")
    ap.add_argument("--name", default="kintcheck_report", help="base name of the report files")
    ap.add_argument("--progress", type=float, default=0, help="print simulated-time progress every N seconds")
    ap.add_argument("--switches", help="override the pack's Kintecus switches")
    ap.add_argument("--no-projection", action="store_true",
                    help="never stop a run early because its projected time exceeds --cap; give it the full --cap")
    ap.add_argument("--keep", action="store_true", help="keep each pack's run folder (default: keep failures only)")
    ap.add_argument("--show", action="store_true", help="show the Kintecus console windows (Windows)")
    ap.add_argument("--runner", help='launch template, e.g. "wine {cmd}" (testing outside Windows)')
    a = ap.parse_args()
    global PROGRESS, NO_PROJECTION
    PROGRESS = a.progress
    NO_PROJECTION = a.no_projection
    job_sw = {}
    if a.jobs:
        with open(a.jobs, newline="", encoding="utf-8-sig") as fj:
            for jr in csv.DictReader(fj):
                if jr.get("id"):
                    job_sw[jr["id"].strip().upper()] = (jr.get("switches") or "").strip()

    exe = os.path.abspath(a.kintecus)
    if not os.path.exists(exe):
        sys.exit(f"Kintecus.exe not found at {exe}  (use --kintecus PATH or set KINTECUS_EXE)")
    packs = sorted(glob.glob(os.path.join(a.packs, "**", "*.kintpack"), recursive=True))
    if a.only:
        want = {x.upper() for x in a.only}
        packs = [p for p in packs if os.path.basename(p)[:-9].upper() in want]
    if a.folder:
        packs = [p for p in packs if a.folder.lower() in os.path.basename(os.path.dirname(p)).lower()]
    if job_sw:
        order = list(job_sw)
        packs = sorted([p for p in packs if os.path.basename(p)[:-9].upper() in job_sw],
                       key=lambda p: order.index(os.path.basename(p)[:-9].upper()))
    if not packs:
        sys.exit("no .kintpack files found")
    os.makedirs(a.work, exist_ok=True)
    rep_csv = os.path.join(a.out, a.name + ".csv")
    fields = ["id", "folder", "check", "match_error_percent", "success", "worst_species", "kintecus_state",
              "kintecus_message", "seconds", "rows_compared", "result_when_built", "error_when_built_percent",
              "switches", "note"]
    fh = open(rep_csv, "w", newline="")
    w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); fh.flush()
    counts, t_all = {}, time.time()
    print(f"kintcheck: {len(packs)} packs, Kintecus = {exe}, success = match error < {a.tol * 100:g}% of peak, "
          f"limits: {'no time limit' if not a.cap else f'{a.cap:g} s per pack'}, "
          f"{'no stall limit' if not a.stall else f'stall {a.stall:g} s'}")
    for i, pk in enumerate(packs, 1):
        mid = os.path.basename(pk)[:-9]
        row = {"id": mid, "folder": os.path.basename(os.path.dirname(pk))}
        try:
            pack = json.load(open(pk, encoding="utf-8"))
            sw = a.switches or job_sw.get(mid.upper()) or pack["settings"].get("switches") or "-ig:mass -INT:1 -obeymaxint"
            row["switches"] = sw
            ref_path = os.path.join(a.refs, mid + ".csv")
            ref = read_reference(ref_path) if os.path.exists(ref_path) else None
            t_end = samp = None
            if ref:
                meta, names, scale, rrows = ref
                t_end, samp = float(meta["t_end"][0]), float(meta["t_end"][2])
                row["result_when_built"] = meta.get("result_here", [""])[0]
                eb = meta.get("result_here", ["", "", ""])[2]
                row["error_when_built_percent"] = "" if not eb else f"{float(eb) * 100:.4g}"
            d = os.path.join(a.work, mid)
            shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
            write_deck(pack, d)
            state, msg, wall = run_kintecus(exe, sw, d, a.runner, t_end or 0, a.stall, a.cap, a.show)
            row.update(kintecus_state=state, kintecus_message=msg[:200], seconds=f"{wall:.1f}")
            conc = os.path.join(d, "CONC.TXT")
            tk, cols = read_conc(conc) if os.path.exists(conc) else ([], {})
            if not ref:
                row["check"] = "NO_REFERENCE"
                row["note"] = "no SBML reference for this pack (the SBML could not be simulated or nothing to compare)"
            elif state in ("fatal", "launch_error"):
                row["check"] = "KINTECUS_FATAL" if state == "fatal" else "LAUNCH_ERROR"
            else:
                err, wsp, note, n = compare(tk, cols, names, scale, rrows, samp)
                complete = bool(tk) and max(tk) >= t_end * (1 - 1e-6)
                row.update(worst_species=wsp or "", rows_compared=n, note=note)
                if err is not None:
                    row["match_error_percent"] = "inf" if not math.isfinite(err) else f"{err * 100:.4g}"
                if not complete:
                    row["check"] = {"stalled": "STALLED", "too_slow": "TOO_SLOW", "timeout": "TIMEOUT"}.get(state, "INCOMPLETE")
                    if err is not None:
                        row["note"] = (note + "; " if note else "") + \
                            f"error over the part that ran (to t={max(tk):g} of {t_end:g})"
                elif err is None:
                    row["check"] = "NOTHING_COMPARED"
                else:
                    row["check"] = "PASS" if err < a.tol else "OFF"
            row["success"] = "YES" if row.get("check") == "PASS" else "NO"
            if row["success"] == "YES" and not a.keep:
                shutil.rmtree(d, ignore_errors=True)
            else:
                trim_big(d)
        except Exception as ex:
            row.update(check="CHECK_ERROR", success="NO", note=str(ex)[:200])
        counts[row["check"]] = counts.get(row["check"], 0) + 1
        w.writerow(row); fh.flush()
        print(f"[{i}/{len(packs)}] {mid:18s} {row['check']:14s} "
              f"{(row.get('match_error_percent') or '-') + '%':>10s}  {row.get('kintecus_message', '')[:60]}", flush=True)
    fh.close()
    lines = [f"kintcheck report  {datetime.datetime.now():%Y-%m-%d %H:%M}",
             f"Kintecus: {exe}", f"packs: {len(packs)}   success threshold: match error < {a.tol * 100:g}% of peak",
             f"total time: {(time.time() - t_all) / 60:.1f} min", ""]
    legend = {"PASS": "Kintecus output matches the SBML reference within the threshold",
              "OFF": "Kintecus ran to the end but differs from the SBML reference by more than the threshold",
              "STALLED": "Kintecus stopped producing output (watchdog --stall)",
              "TOO_SLOW": "Kintecus would take longer than --cap seconds",
              "TIMEOUT": "Kintecus ran longer than --cap seconds",
              "INCOMPLETE": "Kintecus ended before the final time",
              "KINTECUS_FATAL": "Kintecus stopped with a FATAL error (see kintecus_message)",
              "NO_REFERENCE": "no SBML reference shipped for this pack",
              "NOTHING_COMPARED": "no comparable species", "LAUNCH_ERROR": "Kintecus.exe could not be started",
              "CHECK_ERROR": "kintcheck itself failed on this pack (see note)"}
    for k, v in sorted(counts.items(), key=lambda x: -x[1]):
        lines.append(f"  {k:16s} {v:5d}   {legend.get(k, '')}")
    # ---- regressions: packs that passed when built and now fail, or pass with a much larger error
    fails = []
    for row in csv.DictReader(open(rep_csv, newline="")):
        built = (row.get("result_when_built") or "").upper()
        if not built.startswith("PASS"):
            continue
        eb, en = row.get("error_when_built_percent") or "", row.get("match_error_percent") or ""
        if row["success"] != "YES":
            row["regression"] = "FAILED"
        else:
            try:
                b, n = float(eb), float(en)
            except ValueError:
                continue
            if n > max(10 * b, 0.1):
                row["regression"] = "WORSE"
            else:
                continue
        fails.append(row)
    fail_csv = os.path.join(a.out, a.name + "_FAILED.csv")
    with open(fail_csv, "w", newline="") as fh2:
        w2 = csv.DictWriter(fh2, fieldnames=["regression"] + fields, extrasaction="ignore"); w2.writeheader()
        for row in fails:
            w2.writerow(row)
    lines += ["", f"REGRESSIONS (passed when built): {len(fails)}"
              + ("" if fails else "   none - every pack that passed when built still passes"),
              "  FAILED = no longer under the threshold; WORSE = still passes but error > 10x the built error and > 0.1%"]
    if fails:
        lines.append(f"  {'id':18s} {'regression':10s} {'check':14s} {'error % now':>12s} {'when built':>11s}  message")
        for row in fails:
            lines.append(f"  {row['id']:18s} {row['regression']:10s} {row['check']:14s} "
                         f"{(row.get('match_error_percent') or '-'):>12s} {(row.get('error_when_built_percent') or '-'):>11s}  "
                         f"{(row.get('kintecus_message') or row.get('note') or '')[:60]}")
    lines.append(f"  table: {fail_csv}")
    lines += ["", f"Details per pack: {rep_csv}",
              "Failed runs are kept in the work folder for inspection: " + a.work]
    open(os.path.join(a.out, a.name + ".txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
