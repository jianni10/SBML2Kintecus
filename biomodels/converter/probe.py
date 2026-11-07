"""Probe/run harness: build a Kintecus deck from strings, run under Wine on a pty, return CONC as dict of arrays."""
import os, re, signal, shutil, subprocess
import numpy as np

EXE = "/home/claude/kin/Kintecus.exe"
PARM = """#\tParameter Description SpreadSheet\t\t\t\t
# Starting Integration Time\tMaximum Integration Time\tEa UNITS\tConc. Units\tX0\t
0.000001\t{maxint}\tKCAL\tMOLES/LITER\t0\t
#Temperature (K) or Filename\tPressure (Constant ? (Yes/No) )\tVolume Profile\tExt\tX0\t
298\tNOPE\tNO\tNo\t0\t
#       DAYS\tHours\tMinutes\tSeconds\tPicoSeconds\t
0\t0\t0\t{secs}\t0\t
#hv(filename)\tSampling Interval (s)\tPercent(%)\tAccuracy\tX0\t
None1\t{samp}\t0\t{acc}\t0\t
# CSTR/PFR inlet Flow Temperature\tX0\tX0\tX0\tX0\t
0\t0\t0\t0\t0\t
END\t\t\t\t\t
"""
ANSI = re.compile(r"\x1b\[[0-9;?]*[a-zA-Z]")


def _conc_progress(d):
    """(rows, last simulated time) of CONC.TXT so far."""
    try:
        with open(f"{d}/CONC.TXT", "rb") as fh:
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


def kill_strays(d):
    """Wine re-parents Kintecus.exe to init, so killing the pty's process group can leave it running.
    Kill every Kintecus.exe whose working directory is this run directory."""
    real = os.path.realpath(d)
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            if open(f"/proc/{pid}/comm").read().strip() == "Kintecus.exe" and \
                    os.path.realpath(os.readlink(f"/proc/{pid}/cwd")) == real:
                os.kill(int(pid), signal.SIGKILL)
        except (OSError, ProcessLookupError):
            pass


def kintecus(d, switches, timeout=60, t_end=None, stall=15.0, cap=600.0):
    """Run Kintecus.exe -i in directory d on a pseudo-terminal. KRUN.TMP / KDONE.OUT / KERR.OUT report state;
    CONC.TXT growth is the progress signal. Kill on: no new CONC row for `stall` s, a projected finish beyond
    `cap`, or `cap` wall seconds. Returns (console_text, timed_out). `timeout` is kept for callers that pass
    only it (it then acts as the cap when t_end is unknown)."""
    import time
    cmd = ["script", "-qfec", f"/usr/lib/wine/wine64 Kintecus.exe {switches} -i", "/dev/null"]
    out_path = f"{d}/.console.raw"
    with open(out_path, "wb") as fo:
        p = subprocess.Popen(cmd, cwd=d, stdin=subprocess.PIPE, stdout=fo, stderr=subprocess.STDOUT,
                             env={**os.environ, "WINEDEBUG": "-all"}, start_new_session=True)
        try:
            p.stdin.write(b"\r\n" * 400); p.stdin.close()
        except OSError:
            pass
        start = time.time(); last_rows, last_change = 0, start
        limit = cap if t_end else timeout
        killed = None
        while p.poll() is None:
            time.sleep(0.5)
            now = time.time()
            if os.path.exists(f"{d}/KDONE.OUT") or os.path.exists(f"{d}/KERR.OUT"):
                try:
                    p.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    pass
                break
            rows, tsim = _conc_progress(d)
            if rows > last_rows:
                last_rows, last_change = rows, now
            elapsed = now - start
            if now - last_change > stall and elapsed > stall:
                killed = f"stalled (no output row for {stall:.0f} s at t={tsim:g})"
            elif t_end and elapsed > 10 and tsim > 0 and elapsed * t_end / tsim > cap:
                killed = f"too slow (projected {elapsed * t_end / tsim:.0f} s at t={tsim:g})"
            elif elapsed > limit:
                killed = f"exceeded {limit:.0f} s"
            if killed:
                try:
                    os.killpg(p.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                break
        try:
            p.wait(timeout=10)
        except subprocess.TimeoutExpired:
            pass
    kill_strays(d)
    log = ANSI.sub("", open(out_path, "rb").read().decode("latin-1", errors="replace"))
    os.remove(out_path)
    if killed:
        log += f"\n[watchdog] {killed}\n"
    if os.path.exists(f"{d}/KERR.OUT"):            # -i mode reports the fatal number here
        log += "\n[KERR.OUT]\n" + open(f"{d}/KERR.OUT", encoding="latin-1", errors="replace").read()
    open(f"{d}/console.log", "w").write(log)
    return log, bool(killed)


def read_conc(path):
    L = [l.rstrip("\r\n").split("\t") for l in open(path, encoding="latin-1")]
    hdr = [h.strip() for h in L[0] if h.strip()]
    rows = []
    for l in L[1:]:
        if not l or not l[0].strip():
            continue
        try:
            rows.append([float(x) for x in l if x.strip()])
        except ValueError:
            rows.append([float("nan")] * len(hdr))
    D = np.array(rows) if rows else np.zeros((0, len(hdr)))
    return {h: D[:, i] for i, h in enumerate(hdr)}


def run(name, model_lines, species, secs=10, samp=1, maxint=1, acc=1e-8, extra_files=None,
        switches="-INT:1 -obeymaxint", timeout=60):
    d = f"/home/claude/kin/probe/runs/{name}"
    shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    w = lambda f, s: open(f"{d}/{f}", "w", newline="").write(s.replace("\n", "\r\n"))
    w("MODEL.DAT", "# probe\n" + "\n".join(model_lines) + "\nEND\n")
    sp = "# Species\tRes\tInit\tDisp\tExt\tSSD\tConst\tComments\t\n"
    for s in species:   # (name, init, col7)
        sp += f"{s[0]}\t0\t{s[1]}\tY\t0\tNo\t{s[2] if len(s) > 2 else 'No'}\t\t\n"
    w("SPECIES.DAT", sp + "END\n")
    w("PARM.DAT", PARM.format(secs=secs, samp=samp, maxint=maxint, acc=acc))
    for f, s in (extra_files or {}).items():
        w(f, s)
    shutil.copy(EXE, d)
    log, to = kintecus(d, f"-ig:mass {switches}", timeout)
    fat = re.findall(r"FATAL Error #\s*\d+[^\n]*\n[^\n]*\n[^\n]*", log)
    msg = ("TIMEOUT; " if to else "") + (fat[0] if fat else "")
    if not os.path.exists(f"{d}/CONC.TXT"):
        return None, msg or log[-400:]
    return read_conc(f"{d}/CONC.TXT"), msg
