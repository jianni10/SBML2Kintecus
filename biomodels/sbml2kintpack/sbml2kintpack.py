#!/usr/bin/env python3
"""
sbml2kintpack.py -- convert one SBML model into a Kintecus .kintpack (GUI package) plus MODEL.DAT / SPECIES.DAT /
PARM.DAT and any Constant File profiles. Uses the same converter (sbml2kin_batch.py, next to this file) that built
the BioModels collection.

Install once:   pip install python-libsbml sympy            (required)
                pip install libroadrunner numpy              (recommended: exact initial values + automatic duration)

Examples:
    python sbml2kintpack.py BIOMD0000000010.xml
    python sbml2kintpack.py model.xml --t-end 500 --points 1000
    python sbml2kintpack.py model.xml --kintecus "C:\\Kintecus\\Kintecus.exe"     (also checks long USER lines)

Output (in --out, default: next to the SBML file):
    <name>.kintpack              load with the Kintecus GUI ("Load Package...")
    <name>_deck\\                 MODEL.DAT, SPECIES.DAT, PARM.DAT and profile files, ready for Kintecus.exe

Release settings: for a BioModels entry (ID in the file name or the SBML), duration, sampling, switches,
accuracy and deck variant (shifted states etc.) come from release_settings.csv, so the pack matches the verified one (--no-release-settings to skip).

Duration: --t-end if given; else the SED-ML file next to the model (first uniformTimeCourse, unless it is an
auto-generated placeholder); else chosen from the dynamics with libRoadRunner (settled / oscillating / cap);
else 100. Time units are the SBML's, used as Kintecus seconds.

Initial values: with libRoadRunner installed, every species/parameter value at t=0 is taken from it, so initial
assignments and assignment rules come out right. Without it, the values written in the SBML are used and models
that set initial values by rules may start wrong (a warning is printed).

--kintecus: Kintecus V2027 mis-reads some long USER expressions (see KINTECUS_ISSUES.txt #2). With this option
every USER expression longer than 180 characters is evaluated by Kintecus in a probe deck and compared with the
intended value; mis-read lines are re-padded with redundant parentheses until Kintecus reads them correctly.
"""
import argparse, datetime, glob, json, math, os, random, shutil, subprocess, sys, tempfile, time, re
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sbml2kin_batch as K

SWITCHES = "-ig:mass -INT:1 -obeymaxint"


# ------------------------------------------------------------------ duration / sampling
def sedml_timecourse(sbml_path, sedml=None):
    files = [sedml] if sedml else sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(sbml_path)), "*.sedml")))
    for f in files:
        try:
            r = ET.parse(f).getroot()
        except (ET.ParseError, OSError):
            continue
        for e in r.iter():
            if e.tag.endswith("uniformTimeCourse"):
                try:
                    s = float(e.get("outputStartTime", 0)); t = float(e.get("outputEndTime"))
                    n = int(float(e.get("numberOfPoints", e.get("numberOfSteps", 100))))
                    if t > 0 and n > 0:
                        return s, t, n, os.path.basename(f), (e.get("id") or "").startswith("auto_")
                except (TypeError, ValueError):
                    pass
    return None


def rr_load(path):
    import roadrunner
    try:
        rr = roadrunner.RoadRunner(path)
    except RuntimeError as ex:
        # libRoadRunner refuses SBML 'fast' reactions; in BioModels the flag is an export artifact (all reactions
        # fast, or fast synthesis/degradation), so the model is simulated as written with the flag removed
        if "fast" not in str(ex):
            raise
        import libsbml as _ls
        d = _ls.readSBMLFromFile(path)
        for r in d.getModel().getListOfReactions():
            if r.isSetFast() and r.getFast():
                r.setFast(False)
        rr = roadrunner.RoadRunner(_ls.writeSBMLToString(d))
    m = rr.model
    init = {}
    for ids, vals in ((m.getFloatingSpeciesIds(), m.getFloatingSpeciesConcentrations()),
                      (m.getBoundarySpeciesIds(), m.getBoundarySpeciesConcentrations()),
                      (m.getGlobalParameterIds(), m.getGlobalParameterValues()),
                      (m.getCompartmentIds(), m.getCompartmentVolumes())):
        for i, v in zip(ids, vals):
            init[i] = float(v)
    return rr, init


def choose_horizon(rr):
    import numpy as np
    sids = list(rr.model.getFloatingSpeciesIds())
    pick, why = 10.0, "default"
    for T in [10, 30, 100, 300, 1000, 3000, 1e4, 3e4, 1e5]:
        try:
            rr.resetToOrigin()
            rr.timeCourseSelections = ["time"] + [f"[{x}]" for x in sids]
            r = np.asarray(rr.simulate(0, T, 401))
        except Exception:
            break
        Y = r[:, 1:]
        if Y.size == 0 or not np.all(np.isfinite(Y)):
            break
        pick = T
        rng = Y.max(0) - Y.min(0)
        big = rng > 1e-9 * max(np.abs(Y).max(), 1e-300)
        if not big.any():
            why = "no dynamics"; break
        tail = Y[int(0.8 * len(Y)):, big]
        if ((tail.max(0) - tail.min(0)) / rng[big]).max() < 0.01:
            why = "settled"; break
        y = Y[:, int(np.argmax(rng * big))]
        hi = y.min() + 0.5 * (y.max() - y.min())
        if np.sum((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]) & (y[1:-1] > hi)) >= 5:
            why = "oscillating"; break
        why = "still changing at the 1e5 cap"
    rr.resetToOrigin()
    return float(pick), why


def choose_sampling(t_end, n_points, event_times):
    base = t_end / max(n_points, 1)
    samp = base
    pts = [t for t in event_times if 0 < t < t_end]
    if pts:
        g = K.Converter.grid(pts + [t_end])
        if g:
            samp = g / math.ceil(g / base - 1e-9)
    if t_end / samp < 61:
        samp = samp / math.ceil(61 / (t_end / samp))
    return samp


# ------------------------------------------------------------------ Kintecus evaluation check of long USER lines
USER_RE = re.compile(r"(S\[USER;([^:]*):\()(.*)(\)\])")
CHECK_LEN = 180
FUN = {"exp": math.exp, "log": math.log, "log10": math.log10, "sqrt": math.sqrt, "abs": abs, "sin": math.sin,
       "cos": math.cos, "tan": math.tan, "asin": math.asin, "acos": math.acos, "atn": math.atan, "sinh": math.sinh,
       "cosh": math.cosh, "tanh": math.tanh, "floor": math.floor}


def py_value(expr, vals):
    env = dict(FUN); env.update(vals); env["__builtins__"] = {}
    try:
        v = eval(expr.replace("^", "**").lower(), env)
        v = complex(v)
        return v.real if abs(v.imag) <= 1e-12 * max(1.0, abs(v.real)) else float("nan")
    except (ValueError, ZeroDivisionError, OverflowError, TypeError):
        return float("nan")


def run_kintecus(exe, d, runner=None, limit=60):
    for f in ("CONC.TXT", "KDONE.OUT", "KERR.OUT"):
        try:
            os.remove(os.path.join(d, f))
        except OSError:
            pass
    args = [exe] + SWITCHES.split() + ["-i"]
    kw = dict(cwd=d, stdin=subprocess.DEVNULL)
    if runner:
        target = runner.replace("{cmd}", " ".join(args)); kw.update(shell=True, stdout=subprocess.DEVNULL,
                                                                    stderr=subprocess.DEVNULL)
        if os.name != "nt":
            kw["start_new_session"] = True
    else:
        target = args
        if os.name == "nt":
            kw["creationflags"] = subprocess.CREATE_NEW_CONSOLE
            si = subprocess.STARTUPINFO(); si.dwFlags |= subprocess.STARTF_USESHOWWINDOW; si.wShowWindow = 0
            kw["startupinfo"] = si
        else:
            kw.update(stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    p = subprocess.Popen(target, **kw)
    t0 = time.time()
    while time.time() - t0 < limit:
        time.sleep(0.5)
        if any(os.path.exists(os.path.join(d, f)) for f in ("KDONE.OUT", "KERR.OUT")) or p.poll() is not None:
            time.sleep(1.0)
            break
    if p.poll() is None:
        try:
            if runner and os.name != "nt":
                import signal
                os.killpg(p.pid, signal.SIGKILL)
            else:
                p.kill()
        except Exception:
            pass
    if runner and os.name != "nt":                 # Wine: Kintecus.exe outlives its launcher; kill it by folder
        real = os.path.realpath(d)
        for pid in os.listdir("/proc"):
            try:
                if pid.isdigit() and open(f"/proc/{pid}/comm").read().strip() == "Kintecus.exe" and \
                        os.path.realpath(os.readlink(f"/proc/{pid}/cwd")) == real:
                    os.kill(int(pid), 9)
            except OSError:
                pass
    if os.path.exists(os.path.join(d, "KERR.OUT")) or not os.path.exists(os.path.join(d, "CONC.TXT")):
        return None
    lines = [l.rstrip("\r\n").split("\t") for l in open(os.path.join(d, "CONC.TXT"), encoding="latin-1") if l.strip()]
    if len(lines) < 2:                             # Kintecus stopped before writing any rows
        return None
    hdr = [h.strip().upper() for h in lines[0]]
    last = None
    for l in lines[1:]:
        try:
            last = [float(x) for x in l if x.strip()]
        except ValueError:
            pass
    if not last or abs(last[0] - 1.0) > 1e-6:
        return None
    return {h: v for h, v in zip([h for h in hdr if h], last)}


def check_long_lines(conv, exe, runner=None, log=print):
    long_idx = [i for i, (k, rx, cm) in enumerate(conv.lines)
                if USER_RE.search(rx) and len(USER_RE.search(rx).group(3)) > CHECK_LEN]
    if not long_idx:
        return 0, 0, 0
    rnd = random.Random(12345)
    _, species, parm = K.grids(conv, ["evaluation check"], 1.0, 10, 1e-8)
    vals = {}
    for r in species[2:-1]:
        if str(r[5]).upper().startswith("OUTC["):
            continue                              # output-only species (shifted states): not a state, never pinned
        c0 = float(r[2]); v = (abs(c0) if c0 != 0 else 0.01) * rnd.uniform(0.6, 1.4)
        r[2] = r[6] = K.num(v); r[3] = "N"; vals["c" + r[0].lower()] = float(r[2])
    work = tempfile.mkdtemp(prefix="kin_evalcheck_")

    def expected(rx):
        m = USER_RE.search(rx); v = py_value(m.group(3), vals)
        for s in rx[:m.start()].rstrip("+").split("+"):
            s = s.strip()
            if s and s.upper() not in ("NULL", "0 NULL"):
                m_ = re.match(r"^(\d*\.\d+|\d+)(\D.*)$", s)
                if "c" + s.lower() not in vals and m_ and "c" + m_.group(2).lower() in vals:
                    v *= vals["c" + m_.group(2).lower()] ** float(m_.group(1))     # fractional coefficient = order
                else:
                    v *= vals["c" + s.lower()]
        return v

    def evaluate(lines):
        d = tempfile.mkdtemp(dir=work)
        model = [["#", "", "", "evaluation check", ""]] + \
                [[k, "0", "0", rx[:USER_RE.search(rx).end()] + f"==>PRB{j}", cm] for j, (k, rx, cm) in enumerate(lines)] + \
                [["END", "", "", "", ""]]
        sp = species[:-1] + [[f"PRB{j}", "0", "0", "Y", "0", "No", "No", "probe", ""] for j in range(len(lines))] + [species[-1]]
        K.write_deck(d, model, sp, parm, [])
        got = run_kintecus(exe, d, runner)
        shutil.rmtree(d, ignore_errors=True)
        return [None if got is None else got.get(f"PRB{j}") for j in range(len(lines))]

    def good(rx, v):
        e = expected(rx)
        return v is not None and math.isfinite(e) and abs(v - e) <= 2e-5 * abs(e) + 1e-300

    lines = [conv.lines[i] for i in long_idx]
    first = evaluate(lines)
    bad = [i for i, ln, v in zip(long_idx, lines, first) if not good(ln[1], v)]
    if bad and all(v is None for v in first):        # whole probe deck rejected: test one by one
        bad = [i for i in long_idx if not good(conv.lines[i][1], evaluate([conv.lines[i]])[0])]
    fixed = 0
    for i in bad:
        k, rx, cm = conv.lines[i]
        for pad in range(1, 9):
            cand = USER_RE.sub(lambda m: m.group(1) + "(" * pad + m.group(3) + ")" * pad + m.group(4), rx, count=1)
            if good(cand, evaluate([(k, cand, cm)])[0]):
                conv.lines[i] = (k, cand, cm); fixed += 1
                conv.flag("info", f"{cm}: USER expression wrapped in {pad} redundant parentheses (Kintecus mis-read the "
                                  f"unpadded text)")
                break
        else:
            conv.flag("partial", f"{cm}: Kintecus mis-reads this long USER expression in every equivalent form tried")
    shutil.rmtree(work, ignore_errors=True)
    return len(long_idx), len(bad), fixed


# ------------------------------------------------------------------ main
def negative_shifts(rr, path, t_end, n=1000):
    """sid -> shift for every state the SBML run takes below zero by more than 0.1% of its peak: species and
    parameters with rate rules. shift = 1-2-5 step at or above 2 x |min| (so state = value + shift stays > 0)."""
    import libsbml as _ls
    m = _ls.readSBMLFromFile(path).getModel()
    ids = list(rr.model.getFloatingSpeciesIds())
    ids += [r.getVariable() for r in m.getListOfRules() if r.isRate() and r.getVariable() not in ids
            and m.getParameter(r.getVariable()) is not None]
    out = {}
    if not ids:
        return out
    try:
        rr.resetToOrigin()
        rr.timeCourseSelections = ["time"] + [f"[{i}]" if i in rr.model.getFloatingSpeciesIds() else i for i in ids]
        a = __import__('numpy').asarray(rr.simulate(0, t_end, n + 1))
    except Exception:
        return out
    finally:
        try:
            rr.resetToOrigin()
        except Exception:
            pass
    fin = a[:, 1:][__import__('numpy').isfinite(a[:, 1:])]
    gpk = float(__import__('numpy').max(__import__('numpy').abs(fin))) if fin.size else 0.0
    for j, sid in enumerate(ids, 1):
        col = a[:, j]
        if not __import__('numpy').all(__import__('numpy').isfinite(col)):
            continue
        pk, lo = float(__import__('numpy').max(__import__('numpy').abs(col))), float(__import__('numpy').min(col))
        # below 1e-10 of the largest state the negative values are integration noise of a ~zero variable
        if pk > 0 and lo < -1e-3 * pk and lo < -1e-10 * gpk and abs(lo) < 1e15:     # beyond 1e15 the SBML run itself has diverged
            need = 2 * abs(lo)
            e = 10 ** math.floor(math.log10(need))
            out[sid] = next(f * e for f in (1, 2, 5, 10) if f * e >= need)
    return out


def release_setting(path):
    """row of release_settings.csv for a BioModels entry (ID from the file name or the SBML annotation), else None."""
    f = os.path.join(HERE, "release_settings.csv")
    if not os.path.exists(f):
        return None
    mid = model_id(path)
    if not mid:
        return None
    import csv
    for r in csv.DictReader(open(f, newline="")):
        if r["id"] == mid:
            return r
    return None


def model_id(path):
    """BioModels ID of this file: in the file name, else the folder name (biomodels_original_sbml.zip has one folder
    per entry), else the model's own bqmodel:is annotation. Other BIOMD IDs in the file (isDerivedFrom, citations of
    related models) are never used."""
    for name in (os.path.basename(path), os.path.basename(os.path.dirname(os.path.abspath(path)))):
        m = re.fullmatch(r"(BIOMD\d{10})(_url)?(\.xml|\.sbml)?", name, re.I) or re.match(r"(BIOMD\d{10})\b", name)
        if m:
            return m.group(1).upper()
    try:
        import libsbml
        mod = libsbml.readSBMLFromFile(path).getModel()
        for i in range(mod.getNumCVTerms() if mod is not None else 0):
            cv = mod.getCVTerm(i)
            if cv.getQualifierType() == libsbml.MODEL_QUALIFIER and cv.getModelQualifierType() == libsbml.BQM_IS:
                for k in range(cv.getNumResources()):
                    m = re.search(r"BIOMD\d{10}", cv.getResourceURI(k))
                    if m:
                        return m.group(0)
    except Exception:
        pass
    return None


def main():
    ap = argparse.ArgumentParser(description="SBML -> Kintecus .kintpack + DAT files")
    ap.add_argument("sbml", help="SBML file (.xml / .sbml)")
    ap.add_argument("--t-end", type=float, help="simulation end time (SBML time units)")
    ap.add_argument("--points", type=int, help="number of output intervals (default: SED-ML or 1000)")
    ap.add_argument("--sedml", help="SED-ML file to take the duration from (default: any .sedml next to the model)")
    ap.add_argument("--out", help="output folder (default: next to the SBML file)")
    ap.add_argument("--name", help="base name for the outputs (default: SBML file name)")
    ap.add_argument("--switches", help=f"Kintecus switches stored in the pack (default: the release setting for a "
                                       f"BioModels entry, else {SWITCHES})")
    ap.add_argument("--accuracy", type=float, help="PARM.DAT accuracy (default: the release setting for a BioModels "
                                                  "entry, else 1e-8)")
    ap.add_argument("--no-shift", action="store_true", help="do not shift states that go below zero (signed lines only)")
    ap.add_argument("--no-release-settings", action="store_true",
                    help="ignore release_settings.csv (duration, sampling, switches, accuracy of the released packs)")
    ap.add_argument("--no-roadrunner", action="store_true", help="do not use libRoadRunner even if installed")
    ap.add_argument("--kintecus", help="Kintecus.exe: evaluation check of long USER lines (recommended)")
    ap.add_argument("--runner", help='launch template for --kintecus outside Windows, e.g. "wine {cmd}"')
    a = ap.parse_args()

    path = os.path.abspath(a.sbml)
    name = a.name or os.path.splitext(os.path.basename(path))[0]
    out = os.path.abspath(a.out or os.path.dirname(path))
    notes = []
    rr, init = None, {}
    if not a.no_roadrunner:
        try:
            rr, init = rr_load(path)
        except ImportError:
            notes.append("libRoadRunner not installed: initial values taken from the SBML as written "
                         "(models that set them by rules may start wrong)")
        except Exception as ex:
            notes.append(f"libRoadRunner could not load the model ({str(ex)[:120]}); SBML values used")
    # release settings: a BioModels entry converted with the same duration, sampling, switches and accuracy as the
    # verified pack in biomodels_kintpacks_final.zip (release_settings.csv next to this file; also the deck variant: shifted states, unshifted, or legacy)
    rel = None if a.no_release_settings else release_setting(path)
    if rel:
        notes.append(f"{rel['id']}: duration, sampling, switches and accuracy of the released (verified) pack used "
                     f"(--no-release-settings to choose them here)")
        a.switches = a.switches or rel["switches"]
        if a.accuracy is None:
            a.accuracy = float(rel["accuracy"])
    a.switches = a.switches or SWITCHES
    if a.accuracy is None:
        a.accuracy = 1e-8
    # duration
    t_start, n_pts, dur_src = 0.0, a.points or 1000, ""
    tc = sedml_timecourse(path, a.sedml)
    if rel and not a.t_end and not a.points:
        t_end, dur_src = float(rel["t_end"]), f"release setting ({rel['id']})"
    elif a.t_end:
        t_end, dur_src = a.t_end, "given on the command line"
    elif tc and not tc[4]:
        t_start, t_end, n_pts, dur_src = tc[0], tc[1], a.points or tc[2], f"SED-ML {tc[3]}"
    elif rr is not None:
        t_end, why = choose_horizon(rr)
        dur_src = f"chosen from the dynamics ({why})"
    else:
        t_end, dur_src = 100.0, "default (no --t-end, no SED-ML, no libRoadRunner)"
        notes.append("duration defaulted to 100; pass --t-end")
    # states the SBML takes below zero are written shifted (state + shift > 0) with an OUTC[] output of the true value
    variant = rel["variant"] if rel and rel.get("variant") else "plain"
    if rel and variant == "shifted":
        shifts = {k: float(v) for k, v in json.loads(rel["offsets"] or "{}").items()}
    elif rel or rr is None or a.no_shift:
        shifts = {}
    else:
        shifts = negative_shifts(rr, path, t_end)
    if shifts:
        notes.append(f"{len(shifts)} state(s) go below zero in the SBML: written shifted, true value as OUTC[] output")
    if variant == "legacy":                  # deck written as by the converter before 2026-11-01
        K.DIV_MODE, K.ARULE_SWITCH = "eps", False
    conv = K.convert(path, init, t_end, offsets=shifts) if shifts else K.convert(path, init, t_end)
    ev = sorted({e["t"] for e in conv.events} | set(getattr(conv, "switch_times", [])))
    samp = float(rel["sampling"]) if rel and not a.t_end and not a.points else choose_sampling(t_end, n_pts, ev)
    if t_end / samp > 100000:
        samp = choose_sampling(t_end, n_pts, [])
        conv.flag("partial", f"event times do not fit a practical sampling grid; events fire up to one sampling "
                             f"interval ({K.num(samp)}) late")
    if a.kintecus:
        n, b, f = check_long_lines(conv, os.path.abspath(a.kintecus), a.runner)
        if n:
            print(f"evaluation check: {n} long USER lines, {b} mis-read by Kintecus, {f} fixed by re-padding")
    elif any(USER_RE.search(rx) and len(USER_RE.search(rx).group(3)) > CHECK_LEN for _, rx, _ in conv.lines):
        notes.append("long USER expressions present and not checked in Kintecus (use --kintecus PATH)")
    if not conv.lines:
        sys.exit("nothing convertible in this model:\n  " + "\n  ".join(t for s, t in conv.flags if s == "blocker"))
    title = conv.m.getName() or conv.m.getId() or name
    header = [f"{title}  (source file {os.path.basename(path)})", "Converted from SBML by sbml2kintpack.py"]
    model, species, parm = K.grids(conv, header, t_end, int(round(t_end / samp)), a.accuracy)
    deck = os.path.join(out, f"{name}_deck")
    shutil.rmtree(deck, ignore_errors=True); os.makedirs(deck)
    K.write_deck(deck, model, species, parm, conv.ofiles)
    desc = [title, f"SBML file: {os.path.basename(path)}.",
            f"Simulation: 0 to {K.num(t_end)} time units ({dur_src}), sampling {K.num(samp)}.",
            f"Run switches: {a.switches}", "Not verified (converted by sbml2kintpack.py).",
            "SBML time units are used unchanged as Kintecus seconds."]
    if conv.flags or notes:
        desc.append("Conversion notes:")
        desc += [f"  [{s}] {t}" for s, t in conv.flags] + [f"  [note] {t}" for t in notes]
    alias = [f"{r[0]} = {r[5]}" for r in conv.species_rows if r[0] != r[5] and r[5] != "__time__"]
    if alias:
        desc.append("Renamed species (Kintecus name = SBML id): " + "; ".join(alias))
    pack = K.kintpack(model, species, parm, conv.ofiles, a.switches, "\n".join(desc),
                      datetime.datetime.now(datetime.timezone.utc).isoformat())
    kp = os.path.join(out, f"{name}.kintpack")
    json.dump(pack, open(kp, "w"), indent=1)
    blockers = [t for s, t in conv.flags if s == "blocker"]
    print(f"{title}\n  {conv.m.getNumReactions()} SBML reactions -> {len(conv.lines)} Kintecus lines "
          f"({sum('S[USER;' in rx for _, rx, _ in conv.lines)} USER), {len(conv.species_rows)} species")
    print(f"  0..{K.num(t_end)} ({dur_src}), sampling {K.num(samp)}")
    print(f"  wrote {kp}\n  wrote {deck}{os.sep}  (MODEL.DAT, SPECIES.DAT, PARM.DAT"
          f"{', ' + ', '.join(o['filename'] for o in conv.ofiles) if conv.ofiles else ''})")
    print(f"  run:  Kintecus.exe {a.switches}   (in {deck})")
    if blockers:
        print(f"  NOT converted ({len(blockers)}):")
        for t in blockers[:10]:
            print("    - " + t)
    for t in notes:
        print("  note: " + t)


if __name__ == "__main__":
    main()
