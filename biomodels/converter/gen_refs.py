"""gen_refs.py [IDFILE] -- write references/<ID>.csv: the libRoadRunner (SBML) time course of every compared species,
on the pack's own sampling grid (for packs with Time/value-pair files: the rows Kintecus writes,
see kintecus_grid), restricted to the rows the verification compares (t >= SED-ML output start, rows
within half a sampling interval of an event excluded), keyed by Kintecus species name. kintcheck.py compares
Kintecus output against these files: max |Kintecus - reference| / (species peak in the reference)."""
import csv, json, math, os, signal, sys
import numpy as np

B = "/home/claude/kin/batch"
sys.path.insert(0, B)
import run_model as R
import sbml2kin_batch as K
import roadrunner

REF = f"{B}/references"
MAX_ROWS = 1000
os.makedirs(REF, exist_ok=True)


def fmt(v):
    return "nan" if not math.isfinite(v) else f"{v:.9g}"


def pair_times(pack):
    """times of every Time/value-pair Constant File in the pack (first row contains 'Time', pairs until END)."""
    ts = set()
    for o in pack.get("ofiles", []):
        rows = o.get("rows") or []
        if not rows or "TIME" not in str(rows[0][0]).upper():
            continue
        for r in rows[1:]:
            if not r or str(r[0]).strip().upper() == "END":
                break
            ts.add(float(r[0]))
    return sorted(ts)


def kintecus_grid(t_end, samp, evs):
    """Kintecus output rows: 0, samp, 2 samp, ...; at each pair time te not on the grid it writes a row at te and
    restarts the grid from te (te + samp, ...); last row at t_end. Checked against 148 finished pair runs."""
    out, t0, i = [0.0], 0.0, 0
    evs = [e for e in evs if 0 < e < t_end]
    eps = 1e-9 * samp
    while True:
        k = 1
        while True:
            t = t0 + k * samp
            nxt = evs[i] if i < len(evs) else None
            if nxt is not None and abs(t - nxt) <= eps:
                i += 1
            elif nxt is not None and t > nxt:
                out.append(nxt); t0 = nxt; i += 1
                break
            if t >= t_end - eps:
                out.append(t_end)
                return np.array(out)
            out.append(t); k += 1


def make(mid):
    res = json.load(open(f"{R.OUT}/{mid}/result.json"))
    if not res.get("kintpack"):
        return "no pack"
    pack = json.load(open(f"{R.PACKS}/{res['kintpack']}"))
    path = f"{R.ROOT}/{mid}/{R.MASTERS[mid]}"
    rr, init = R.rr_load(path)
    t_end, samp, t_start = res["t_end"], res["sampling"], res.get("t_start", 0.0) or 0.0

    def _alarm(*_):
        raise TimeoutError("conversion over 180 s")
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(180)
    try:
        conv = K.convert(path, init, t_end)
    finally:
        signal.alarm(0)
    pack_names = {r[0].upper() for r in pack["grids"]["species"]["rows"] if r[0] and not r[0].startswith("#")}
    names = [(r[0], r[5]) for r in conv.species_rows if r[2] == "Y"]
    if not {n.upper() for n, _ in names} <= pack_names:
        # pack built by an older converter: take Kintecus name <-> SBML id from the pack ("label [sid]" comments)
        import re as _re
        names = []
        for r in pack["grids"]["species"]["rows"]:
            m = _re.search(r"\[([^\[\]]+)\]\s*$", r[7] if len(r) > 7 else "")
            if r[0] and not r[0].startswith("#") and r[0] != "END" and r[3].upper().startswith("Y") and m:
                names.append((r[0], m.group(1)))
        if not names:
            return "species names differ from the pack and the pack has no id comments"
    ev = sorted({e["t"] for e in conv.events} | set(getattr(conv, "switch_times", [])))
    times = np.round(np.linspace(0, t_end, int(round(t_end / samp)) + 1), 12)
    ptimes = pair_times(pack)
    if ptimes:                                        # Kintecus 2026-11-05: output grid restarts at each pair time
        times = np.round(kintecus_grid(t_end, samp, ptimes), 12)
    rr.resetToOrigin()
    sids = set(rr.model.getFloatingSpeciesIds()) | set(rr.model.getBoundarySpeciesIds())
    gps = set(rr.model.getGlobalParameterIds())
    sel = [(n, s) for n, s in names if s in sids or s in gps]
    if not sel:
        return "nothing to compare"
    rr.timeCourseSelections = ["time"] + [f"[{s}]" if s in sids else s for _, s in sel]
    DEFAULT_TOL = os.environ.get("REF_DEFAULT_TOL") == "1"   # fallback for models CVODE cannot do at 1e-10
    try:
        if not DEFAULT_TOL:
            rr.integrator.relative_tolerance = 1e-10
            rr.integrator.absolute_tolerance = 1e-16
            rr.integrator.maximum_num_steps = 200000
    except Exception:
        pass
    ref = np.asarray(rr.simulate(times=list(times)))
    keep = times >= t_start - 1e-12
    for te in ev:
        keep &= np.abs(times - te) > 0.5 * samp
    if not keep.any():
        return "no comparable rows"
    scales = [float(np.nanmax(np.abs(ref[keep, j]))) for j in range(1, len(sel) + 1)]
    gmax = max(scales)
    cols = [j for j, sc in enumerate(scales, 1) if sc > 1e-12 * gmax and sc > 0]
    if not cols:
        return "all compared species stay at zero"
    rows = np.flatnonzero(keep)
    pick = set(np.unique(np.round(np.linspace(0, len(rows) - 1, min(MAX_ROWS, len(rows)))).astype(int)))
    for j in cols:                                    # always keep each species' extreme points
        v = ref[rows, j]
        pick |= {int(np.nanargmax(v)), int(np.nanargmin(v))}
    rows = rows[sorted(pick)]
    st = res["status"]
    with open(f"{REF}/{mid}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["#kintcheck-reference", "v1", mid])
        w.writerow(["#t_end", fmt(t_end), "sampling", fmt(samp), "switches", res.get("switches") or ""])
        w.writerow(["#result_here", st, "max_rel_error", "" if res.get("error") is None else fmt(res["error"]),
                    "deck_check", fmt((res.get("deck_check") or {}).get("err", float("nan")))])
        w.writerow(["#source", f"libRoadRunner {roadrunner.__version__}",
                    *(["default tolerances (CVODE fails at rtol 1e-10)", ""] if os.environ.get("REF_DEFAULT_TOL") == "1"
                      else ["rtol 1e-10", "atol 1e-16"]),
                    f"{len(rows)} of {int(keep.sum())} compared rows"])
        w.writerow(["#scale"] + [fmt(scales[j - 1]) for j in cols])
        w.writerow(["Time"] + [sel[j - 1][0] for j in cols])
        for i in rows:
            w.writerow([fmt(times[i])] + [fmt(ref[i, j]) for j in cols])
    return "ok"


if __name__ == "__main__":
    ids = open(sys.argv[1]).read().split() if len(sys.argv) > 1 else \
        sorted(d for d in os.listdir(R.OUT) if os.path.exists(f"{R.OUT}/{d}/result.json"))
    skip = set(open(sys.argv[2]).read().split()) if len(sys.argv) > 2 else set()
    log = open(f"{B}/gen_refs.log", "a")
    for mid in ids:
        if mid in skip:
            continue
        try:
            msg = make(mid)
        except Exception as ex:
            msg = "error: " + str(ex)[:150]
        print(mid, msg, file=log, flush=True)
