"""
run_model.py BIOMD_ID  -- convert one curated BioModels entry, run Kintecus, verify against libRoadRunner, write .kintpack.
Writes out/<ID>/ (deck, CONC.TXT, console logs) and out/<ID>/result.json, and packs/<ID>.kintpack.
"""
import json, math, os, re, shutil, sys, time, glob, datetime, traceback
import xml.etree.ElementTree as ET
from fractions import Fraction
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/home/claude/kin/probe")
import sbml2kin_batch as K
import deck_sim
import signal
from probe import kintecus, read_conc

ROOT = "/home/claude/tbm/final"
OUT = os.environ.get("KIN_OUT", "/home/claude/kin/batch/out")
PACKS = os.environ.get("KIN_PACKS", "/home/claude/kin/batch/packs")
MASTERS = json.load(open("/home/claude/kin/all_masters.json"))
KIN_TIMEOUT = int(os.environ.get("KIN_TIMEOUT", "120"))
PASS, APPROX = 4e-2, 4e-2   # PASS = Kintecus under 4% of the SBML (James, 2026-10-30)
DECK_TOL = 1e-3            # deck verification (independent SciPy run) stays at 0.1%


def sedml_timecourse(mid):
    """First uniformTimeCourse in the entry's SED-ML: (start, end, points)."""
    for f in sorted(glob.glob(f"{ROOT}/{mid}/*.sedml")):
        try:
            r = ET.parse(f).getroot()
        except ET.ParseError:
            continue
        for e in r.iter():
            if e.tag.endswith("uniformTimeCourse"):
                try:
                    s = float(e.get("outputStartTime", 0)); t = float(e.get("outputEndTime"))
                    n = int(float(e.get("numberOfPoints", e.get("numberOfSteps", 100))))
                    if t > 0 and n > 0:
                        placeholder = (e.get("id") or "").startswith("auto_")
                        return s, t, n, os.path.basename(f) + (" (auto-generated placeholder)" if placeholder else "")
                except (TypeError, ValueError):
                    pass
    return 0.0, 100.0, 100, None


def choose_horizon(rr, sids):
    """For entries without a curated time course: smallest T in a 1-3-10 ladder at which the floating species have
    settled (change over the last 20% of [0,T] < 1% of their range) or at least 5 oscillation peaks have occurred."""
    ladder = [10, 30, 100, 300, 1000, 3000, 1e4, 3e4, 1e5]
    pick, why = 10.0, "default"
    for T in ladder:
        try:
            rr.resetToOrigin()
            fs = set(rr.model.getFloatingSpeciesIds()) | set(rr.model.getBoundarySpeciesIds())
            rr.timeCourseSelections = ["time"] + [f"[{x}]" if x in fs else x for x in sids]
            r = np.asarray(rr.simulate(0, T, 401))
        except Exception:
            break
        Y = r[:, 1:]
        if Y.size == 0:
            why = "no floating species"; break
        if not np.all(np.isfinite(Y)):
            break
        pick = T
        rng = Y.max(0) - Y.min(0)
        big = rng > 1e-9 * max(np.abs(Y).max(), 1e-300)
        if not big.any():
            why = "no dynamics"; break
        tail = Y[int(0.8 * len(Y)):, big]
        rel = (tail.max(0) - tail.min(0)) / rng[big]
        if rel.max() < 0.01:
            why = "settled"; break
        j = int(np.argmax(rng * big)); y = Y[:, j]
        hi = y.min() + 0.5 * (y.max() - y.min())
        peaks = np.sum((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]) & (y[1:-1] > hi))
        if peaks >= 5:
            why = "oscillating"; break
        why = "still changing at cap" if T == ladder[-1] else why
    rr.resetToOrigin()
    return float(pick), why


FAST_STRIPPED = set()


def rr_load(path):
    """libRoadRunner refuses SBML 'fast' reactions. In the BioModels entries that use them the flag cannot mean
    rapid equilibrium (all reactions fast, or fast synthesis/degradation: speeding them up makes the model diverge),
    so the reference is the model as written with the flag removed - which is also how the converter writes it."""
    import roadrunner
    try:
        rr = roadrunner.RoadRunner(path)
    except RuntimeError as ex:
        if "fast" not in str(ex):
            raise
        import libsbml as _ls
        d = _ls.readSBMLFromFile(path)
        for r in d.getModel().getListOfReactions():
            if r.isSetFast() and r.getFast():
                r.setFast(False)
        rr = roadrunner.RoadRunner(_ls.writeSBMLToString(d))
        FAST_STRIPPED.add(path)
    m = rr.model
    init = {}
    for ids, vals in ((m.getFloatingSpeciesIds(), m.getFloatingSpeciesConcentrations()),
                      (m.getBoundarySpeciesIds(), m.getBoundarySpeciesConcentrations()),
                      (m.getGlobalParameterIds(), m.getGlobalParameterValues()),
                      (m.getCompartmentIds(), m.getCompartmentVolumes())):
        for i, v in zip(ids, vals):
            init[i] = float(v)
    return rr, init


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
        a = np.asarray(rr.simulate(0, t_end, n + 1))
    except Exception:
        return out
    finally:
        try:
            rr.resetToOrigin()
        except Exception:
            pass
    fin = a[:, 1:][np.isfinite(a[:, 1:])]
    gpk = float(np.max(np.abs(fin))) if fin.size else 0.0
    for j, sid in enumerate(ids, 1):
        col = a[:, j]
        if not np.all(np.isfinite(col)):
            continue
        pk, lo = float(np.max(np.abs(col))), float(np.min(col))
        # below 1e-10 of the largest state the negative values are integration noise of a ~zero variable
        if pk > 0 and lo < -1e-3 * pk and lo < -1e-10 * gpk and abs(lo) < 1e15:     # beyond 1e15 the SBML run itself has diverged
            need = 2 * abs(lo)
            e = 10 ** math.floor(math.log10(need))
            out[sid] = next(f * e for f in (1, 2, 5, 10) if f * e >= need)
    return out


def choose_sampling(t_end, n_points, event_times):
    base = t_end / max(n_points, 1)
    samp = base
    pts = [t for t in event_times if 0 < t < t_end]
    if pts:
        g = K.Converter.grid(pts + [t_end])
        if g and math.ceil(g / base - 1e-9) >= 1:
            samp = g / math.ceil(g / base - 1e-9)
        elif g:
            samp = g          # event times share only a tiny common step (e.g. 81.324245310488 and 55.720620288104)
    if t_end / samp < 61:
        samp = samp / math.ceil(61 / (t_end / samp))
    return samp


def compare(rr, conv, kin, t_start, event_times, samp):
    names = {r[0].upper(): r[5] for r in conv.species_rows if r[2] == "Y"}
    cols = {k.upper(): v for k, v in kin.items()}
    t = cols["TIME(S)"]
    t, idx = np.unique(np.round(t, 12), return_index=True)
    sel = [(sid, k) for k, sid in names.items() if k in cols]
    rr.resetToOrigin()
    species_ids = set(rr.model.getFloatingSpeciesIds()) | set(rr.model.getBoundarySpeciesIds())
    sel = [(sid, k) for sid, k in sel if sid in species_ids or sid in set(rr.model.getGlobalParameterIds())]
    rr.timeCourseSelections = ["time"] + [f"[{sid}]" if sid in species_ids else sid for sid, _ in sel]
    try:
        rr.integrator.relative_tolerance = 1e-10
        rr.integrator.absolute_tolerance = 1e-16
        rr.integrator.maximum_num_steps = 200000
    except Exception:
        pass
    times = t if t[0] == 0 else np.concatenate([[0.0], t])
    try:
        ref = rr.simulate(times=list(times))
    except Exception as ex:
        if "CVODE" not in str(ex):
            raise
        # CVODE cannot run some models at rtol 1e-10 / atol 1e-16 (CV_TOO_MUCH_WORK, CV_CONV_FAILURE): use the
        # libRoadRunner defaults, as gen_refs.py does for the same models (reference header says so)
        rr.resetToOrigin()
        rr.integrator.relative_tolerance = 1e-6
        rr.integrator.absolute_tolerance = 1e-12
        compare.default_tol = True
        ref = rr.simulate(times=list(times))
    ref = np.asarray(ref)[-len(t):, :]
    keep = t >= t_start - 1e-12
    for te in event_times:
        keep &= np.abs(t - te) > 0.5 * samp
    scales = {}
    for j, (sid, k) in enumerate(sel, 1):
        scales[sid] = float(np.nanmax(np.abs(ref[keep, j]))) if keep.any() else 0.0
    gmax = max(scales.values()) if scales else 0.0
    per, worst, worst_sp = {}, 0.0, None
    for j, (sid, k) in enumerate(sel, 1):
        sc = scales[sid]
        if sc <= 1e-12 * gmax or sc == 0:
            continue
        a = cols[k][idx][keep]
        b = ref[keep, j]
        if np.any(~np.isfinite(a)):
            e = float("inf")
        else:
            e = float(np.max(np.abs(a - b)) / sc)
        per[sid] = e
        if e > worst or worst_sp is None:
            worst, worst_sp = e, sid
    return worst, worst_sp, per, len(per), int(keep.sum())


USER_RE = re.compile(r"(S\[USER;([^:]*):\()(.*)(\)\])")
CHECK_LEN = 180          # USER expressions longer than this get the Kintecus evaluation check


def _py_value(expr, vals):
    import math as _m
    env = dict(deck_sim.FUN); env.update(vals); env["__builtins__"] = {}
    with np.errstate(all="ignore"):
        v = eval(expr.replace("^", "**").lower(), env)
    v = complex(v)
    return v.real if abs(v.imag) <= 1e-12 * max(1.0, abs(v.real)) else float("nan")


def accept_long_lines(conv, d):
    """Kintecus V2027 mis-reads some long USER expressions depending only on how the characters fall: FATAL #999,
    or a silently wrong value. Every long USER line is evaluated by Kintecus in a probe deck (all species pinned at
    random values, the line's products replaced by a probe species, so probe(t=1) = rate) and compared with the
    same text evaluated in Python. Lines Kintecus rejects or gets wrong are re-padded with redundant parentheses
    until Kintecus evaluates them correctly."""
    long_idx = []
    for i, (k, rx, cm) in enumerate(conv.lines):
        m = USER_RE.search(rx)
        if m and len(m.group(3)) > CHECK_LEN:
            long_idx.append(i)
    if not long_idx:
        return 0, 0, 0
    rng = np.random.default_rng(12345)
    _, species, parm = K.grids(conv, ["evaluation check"], 1.0, 10, 1e-8)
    vals = {}
    for r in species[2:-1]:
        if str(r[5]).upper().startswith("OUTC["):
            continue                              # output-only species: not a state
        c0 = float(r[2])
        v = (abs(c0) if c0 != 0 else 0.01) * rng.uniform(0.6, 1.4)
        r[6] = K.num(v); r[2] = K.num(v); r[3] = "N"
        vals["c" + r[0].lower()] = float(r[2])
    n = [0]

    def expected(rx):
        m = USER_RE.search(rx)
        lhs = rx[:m.start()].rstrip("+")
        v = _py_value(m.group(3), vals)
        for s in lhs.split("+"):
            s = s.strip()
            if s and s.upper() not in ("NULL", "0 NULL"):
                m_ = re.match(r"^(\d*\.\d+|\d+)(\D.*)$", s)
                if "c" + s.lower() not in vals and m_ and "c" + m_.group(2).lower() in vals:
                    v *= vals["c" + m_.group(2).lower()] ** float(m_.group(1))     # fractional coefficient = order
                else:
                    v *= vals["c" + s.lower()]
        return v

    def probe_rx(rx, j):
        m = USER_RE.search(rx)
        return rx[:m.end()] + f"==>PRB{j}"

    def evaluate(lines):
        """-> list of (accepted, value) per line; one Kintecus run."""
        n[0] += 1
        ad = f"{d}/evalcheck_{n[0]}"; os.makedirs(ad)
        model = [["#", "", "", "evaluation check", ""]] + \
                [[k, "0", "0", probe_rx(rx, j), cm] for j, (k, rx, cm) in enumerate(lines)] + [["END", "", "", "", ""]]
        sp_rows = species[:-1] + [[f"PRB{j}", "0", "0", "Y", "0", "No", "No", "probe", ""] for j in range(len(lines))] \
            + [species[-1]]
        K.write_deck(ad, model, sp_rows, parm, [])
        shutil.copy(K_EXE, ad)
        kintecus(ad, "-ig:mass -INT:1 -obeymaxint", timeout=60, stall=30)
        out = [(False, None)] * len(lines)
        if not os.path.exists(f"{ad}/KERR.OUT") and os.path.exists(f"{ad}/CONC.TXT"):
            c = read_conc(f"{ad}/CONC.TXT")
            t = c.get("Time(s)")
            if t is not None and len(t) and abs(t[-1] - 1.0) < 1e-6:
                out = [(True, float(c[f"PRB{j}"][-1])) if f"PRB{j}" in c else (False, None) for j in range(len(lines))]
        shutil.rmtree(ad, ignore_errors=True)
        return out

    def good(rx, got):
        acc, v = got
        e = expected(rx)
        if not acc or v is None or not math.isfinite(e):
            return False, e
        return abs(v - e) <= 2e-5 * abs(e) + 1e-300, e

    lines = [conv.lines[i] for i in long_idx]
    first = evaluate(lines)
    bad = [i for i, ln, got in zip(long_idx, lines, first) if not good(ln[1], got)[0]]
    if bad and not any(a for a, _ in first):          # whole probe deck rejected: test lines one by one
        bad = [i for i in long_idx if not good(conv.lines[i][1], evaluate([conv.lines[i]])[0])[0]]
    fixed = 0
    for i in bad:
        k, rx, cm = conv.lines[i]
        for pad in range(1, 9):
            cand = USER_RE.sub(lambda m: m.group(1) + "(" * pad + m.group(3) + ")" * pad + m.group(4), rx, count=1)
            ok, e = good(cand, evaluate([(k, cand, cm)])[0])
            if ok:
                conv.lines[i] = (k, cand, cm); fixed += 1
                conv.flag("info", f"{cm}: USER expression wrapped in {pad} redundant parentheses (Kintecus V2027 "
                                  f"mis-read the unpadded text)")
                break
        else:
            conv.flag("partial", f"{cm}: Kintecus V2027 mis-reads this long USER expression (FATAL #999 or a wrong "
                                 f"value) in every equivalent form tried")
    return len(long_idx), len(bad), fixed

def main(mid):
    t0 = time.time()
    path = f"{ROOT}/{mid}/{MASTERS[mid]}"
    d = f"{OUT}/{mid}"
    prev = {}                                  # previous deck + its independent check, reused if the deck is unchanged
    try:
        _r = json.load(open(f"{d}/result.json"))
        _att = ([_r["attempts"][_r["best_run"] - 1]] if _r.get("best_run") and _r["best_run"] <= len(_r["attempts"]) else
                [a for a in _r["attempts"] if a["switches"] == _r.get("switches")])
        _rd = f"{d}/run_{_r['attempts'].index(_att[0]) + 1}"
        if (_r.get("deck_check") or {}).get("err") is not None:
            prev = {"check": _r["deck_check"], "files": {fn: open(f"{_rd}/{fn}", "rb").read()
                                                       for fn in os.listdir(_rd) if fn.upper().endswith((".DAT", ".TXT"))
                                                       and fn.upper() not in ("CONC.TXT", "VIEW.TXT", "ADDSPEC.TXT")}}
    except Exception:
        prev = {}
    shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    res = {"id": mid, "file": MASTERS[mid], "status": "ERROR", "flags": [], "attempts": []}
    compare.default_tol = False
    try:
        s_start, t_end, n_pts, sedml = sedml_timecourse(mid)
        res.update(sedml=sedml, t_start=s_start, t_end=t_end, n_points=n_pts)
        rr = None
        try:
            rr, init = rr_load(path)
            if path in FAST_STRIPPED:
                res["flags"].append(["info", "SBML 'fast' reactions: libRoadRunner cannot simulate them, so the reference "
                                             "is the model with the fast flag removed (ordinary kinetics, as the pack "
                                             "writes them); here the flag cannot mean rapid equilibrium (all reactions "
                                             "fast, or fast synthesis/degradation)"])
        except Exception as ex:
            init = {}
            res["flags"].append(["partial", f"libRoadRunner could not load the model ({str(ex)[:150]}); unverified"])
        if rr is not None and (sedml is None or "placeholder" in sedml):
            fl = list(rr.model.getFloatingSpeciesIds())
            try:
                import libsbml as _ls
                _doc = _ls.readSBMLFromFile(path)
                fl += [r.getVariable() for r in _doc.getModel().getListOfRules()
                       if r.isRate() and r.getVariable() not in fl]
            except Exception:
                pass
            t_end, why = choose_horizon(rr, fl)
            try:
                import libsbml as _ls
                _doc2 = _ls.readSBMLFromFile(path)
                _m = _doc2.getModel()
                _tk = []
                for _e in _m.getListOfEvents():
                    _a = _e.getTrigger().getMath() if _e.isSetTrigger() else None
                    if _a is not None and _a.getNumChildren() == 2:
                        for _c in (_a.getChild(0), _a.getChild(1)):
                            if _c.isNumber():
                                _tk.append(float(_c.getValue()))
                if _tk and max(_tk) * 2 > t_end:
                    t_end = min(1e5, 10 ** math.ceil(math.log10(max(_tk) * 2)))
                    why += f"; extended past the last time event ({K.num(max(_tk))})"
            except Exception:
                pass
            s_start, n_pts = 0.0, 1000
            res["t_end"], res["n_points"], res["horizon_rule"] = t_end, n_pts, why
            res["flags"].append(["info", f"no curated SED-ML time course; duration {K.num(t_end)} chosen from the "
                                         f"dynamics ({why})"])
        def _conv_alarm(*_):
            raise K.Unsupported("conversion took longer than 120 s (expression too large)")
        signal.signal(signal.SIGALRM, _conv_alarm); signal.alarm(120)
        shifts = {}
        if rr is not None and os.environ.get("NEG_SHIFT", "1") == "1":
            shifts = negative_shifts(rr, path, t_end)
            if shifts:
                res["shifted"] = shifts
        try:
            conv = K.convert(path, init, t_end, offsets=shifts)
        finally:
            signal.alarm(0)
        flag_at = len(res["flags"])
        res["flags"] += [list(f) for f in conv.flags]
        name = conv.m.getName() or conv.m.getId() or mid
        res["name"] = name
        res["n_reactions_sbml"] = conv.m.getNumReactions()
        res["n_lines"] = len(conv.lines)
        res["n_user_lines"] = sum("S[USER;" in rx for _, rx, _ in conv.lines)
        res["n_species"] = len(conv.species_rows)
        ev_times = sorted({e["t"] for e in conv.events} | set(getattr(conv, "switch_times", [])))
        samp = choose_sampling(t_end, n_pts, ev_times)
        if t_end / samp > 100000:
            samp = choose_sampling(t_end, n_pts, [])
            conv.flag("partial", "event times do not fit a practical sampling grid; events fire up to one sampling "
                                 f"interval ({K.num(samp)}) late")
        res["sampling"] = samp
        header = [f"{name}  ({mid}, curated BioModels; source file {MASTERS[mid]})",
                  "Converted from SBML by sbml2kin_batch.py; verified against libRoadRunner (see .kintpack description)"]
        blockers = [f for f in conv.flags if f[0] == "blocker"]
        if not conv.lines:
            raise K.Unsupported("nothing convertible")
        nlong, nbad, nfixed = accept_long_lines(conv, d)
        res["long_user_lines"], res["long_user_lines_misread"], res["long_user_lines_repadded"] = nlong, nbad, nfixed
        best = None
        tries = [(1e-8, "-ig:mass -INT:1 -obeymaxint"), (1e-10, os.environ.get("KIN_SW2", "-ig:mass -INT:5 -obeymaxint"))]
        if any("(consumes " in cm or "written once as signed" in cm for _, _, cm in conv.lines) or \
                any(rx.startswith("source==>") and str(k).startswith("-") for k, rx, _ in conv.lines):
            # signed lines: the SBML lets these variables go below zero; DASPK without its non-negativity
            # constraint (-INT:5:10:0) integrates through zero where the other integrators stall
            tries.append((1e-10, "-ig:mass -INT:5:10:0 -obeymaxint"))
        conv_main, legacy_done = conv, False
        tries = [(a, b, conv) for a, b in tries]
        while tries:
            acc, sw, conv = tries.pop(0)
            model, species, parm = K.grids(conv, header, t_end, int(round(t_end / samp)), acc)
            ad = f"{d}/run_{len(res['attempts']) + 1}"; os.makedirs(ad)
            K.write_deck(ad, model, species, parm, conv.ofiles)
            shutil.copy(K_EXE, ad)
            ts = time.time()
            multi_user = res["n_user_lines"] > 1
            log, to = kintecus(ad, sw, t_end=t_end, stall=float(os.environ.get("KIN_STALL", "15")), cap=int(os.environ.get("KIN_CAP", "600")))
            wd = re.findall(r"\[watchdog\] (.*)", log)
            att = {"switches": sw, "accuracy": acc, "wall": round(time.time() - ts, 2), "timeout": to}
            if getattr(conv, "legacy", False):
                att["converter"] = "legacy deck (before 2026-11-01)"
            elif getattr(conv, "unshifted", False):
                att["converter"] = "deck without shifted states"
            if wd:
                att["watchdog"] = wd[0]
            fat = re.findall(r"FATAL Error #\s*(\d+)", log)
            if fat:
                att["fatal"] = fat[0]
                att["fatal_text"] = " ".join(log[log.find("FATAL Error"):].split("\n")[1:4]).strip()[:200]
            sw_err = re.findall(r"SOFTWARE ERROR AT\s*(\S+)", log)
            if sw_err and not fat:
                att["fatal"] = f"SOFTWARE ERROR AT {sw_err[0].rstrip('!')}"
                att["fatal_text"] = "Kintecus internal error; it stopped before integrating"
            os.remove(f"{ad}/Kintecus.exe")
            if os.path.exists(f"{ad}/CONC.TXT") and not fat:
                kin = read_conc(f"{ad}/CONC.TXT")
                tk = kin.get("Time(s)", np.array([]))
                att["t_reached"] = float(tk[-1]) if len(tk) else 0.0
                complete = len(tk) and tk[-1] >= t_end * (1 - 1e-6)
                if complete and re.search(r"\bNaN\b|Infinity", open(f"{ad}/CONC.TXT", errors="replace").read()):
                    complete = False      # Kintecus reached t_end but the output went NaN (only Warning #6 reported)
                    att["nan"] = "NaN in CONC.TXT (Kintecus ran to the end with only Warning #6)"
                if rr is not None and complete:
                    try:
                        err, wsp, per, ncmp, nrows = compare(rr, conv, kin, s_start, ev_times, samp)
                        att.update(err=err, worst_species=wsp, n_compared=ncmp, rows_compared=nrows,
                                   per_species={k: round(v, 8) for k, v in sorted(per.items(), key=lambda x: -x[1])[:10]})
                    except Exception as ex:
                        att["compare_error"] = str(ex)[:200]
                att["complete"] = bool(complete)
            res["attempts"].append(att)
            score = att.get("err", float("inf")) if att.get("complete") else float("inf")
            if att.get("complete") and not math.isfinite(score) and "err" not in att:
                score = 1e300     # ran to the end but no reference to score against: still beats an incomplete run
            # an alternative deck (unshifted / before 11-01) replaces the current one only if it passes
            alt_ok = conv is conv_main or (math.isfinite(score) and score < PASS) or best is None
            if alt_ok and (best is None or score < best[0]):
                best = (score, att, model, species, parm, sw, conv)
            if score <= DECK_TOL:     # good enough: stop trying integrators (PASS itself is < 4%)
                break
            if not tries and not legacy_done and os.environ.get("LEGACY_FALLBACK", "1") == "1" and \
                    not (best and math.isfinite(best[0]) and best[0] < PASS):
                # nothing passed with the current converter: also try the deck the converter wrote before
                # 2026-11-01 ([A]+1e-30 divided forms, abs() time steps in assignment rules) and keep the better one
                legacy_done = True
                if shifts:
                    # the shifted deck did not pass: also try the same converter without shifts
                    try:
                        signal.signal(signal.SIGALRM, _conv_alarm); signal.alarm(120)
                        try:
                            conv_n = K.convert(path, init, t_end)
                        finally:
                            signal.alarm(0)
                        if conv_n.lines:
                            accept_long_lines(conv_n, d)
                            conv_n.unshifted = True
                            tries += [(1e-8, "-ig:mass -INT:1 -obeymaxint", conv_n),
                                      (1e-10, "-ig:mass -INT:5 -obeymaxint", conv_n),
                                      (1e-10, "-ig:mass -INT:5:10:0 -obeymaxint", conv_n)]
                    except Exception as ex:
                        res["flags"].append(["info", f"unshifted deck not built ({str(ex)[:100]})"])
                try:
                    old_mode, old_sw = K.DIV_MODE, K.ARULE_SWITCH
                    K.DIV_MODE, K.ARULE_SWITCH = "eps", False
                    signal.signal(signal.SIGALRM, _conv_alarm); signal.alarm(120)
                    try:
                        conv_l = K.convert(path, init, t_end)
                    finally:
                        signal.alarm(0)
                        K.DIV_MODE, K.ARULE_SWITCH = old_mode, old_sw
                    if [l[1] for l in conv_l.lines] != [l[1] for l in conv_main.lines] and conv_l.lines:
                        accept_long_lines(conv_l, d)
                        conv_l.legacy = True
                        tries += [(1e-8, "-ig:mass -INT:1 -obeymaxint", conv_l), (1e-10, "-ig:mass -INT:5 -obeymaxint", conv_l)]
                except Exception as ex:
                    res["flags"].append(["info", f"legacy deck not built ({str(ex)[:100]})"])

        score, att, model, species, parm, sw, conv = best
        res["best_run"] = res["attempts"].index(att) + 1      # run_<n> that the pack and the score come from
        if conv is not conv_main:                 # the legacy deck won: its flags and counts describe the pack
            res.pop("shifted", None)
            res["flags"][flag_at:flag_at + len(conv_main.flags)] = [list(f) for f in conv.flags]
            res["flags"].append(["info", "the converter's current deck did worse here; this pack is the deck written "
                                         "before 2026-11-01 ([A]+1e-30 divided forms, abs() time steps)"
                                 if getattr(conv, "legacy", False) else
                                 "the deck with shifted states did worse here; this pack is the deck without them"])
            res["n_lines"] = len(conv.lines)
            res["n_user_lines"] = sum("S[USER;" in rx for _, rx, _ in conv.lines)
            res["n_species"] = len(conv.species_rows)
            blockers = [f for f in conv.flags if f[0] == "blocker"]
        if score >= 1e300:
            score = float("inf")      # the "complete, unscored" placeholder
        if getattr(compare, "default_tol", False):
            res["flags"].append(["info", "libRoadRunner could not run this model at rtol 1e-10 / atol 1e-16 (CVODE "
                                         "error); the SBML reference uses its default tolerances (1e-6 / 1e-12)"])
        deck_err = None
        if rr is not None and not (math.isfinite(score) and score <= DECK_TOL):
            run_dir = f"{d}/run_{res['attempts'].index(att) + 1}"
            def _alarm(*_):
                raise TimeoutError(f"deck simulation over {os.environ.get('KIN_DECK_ALARM', '180')} s")
            same = bool(prev) and all(os.path.exists(f"{run_dir}/{fn}") and open(f"{run_dir}/{fn}", "rb").read() == b
                                      for fn, b in prev["files"].items()) and \
                {fn for fn in os.listdir(run_dir) if fn.upper().endswith(".DAT")} <= set(prev["files"])
            signal.signal(signal.SIGALRM, _alarm); signal.alarm(int(os.environ.get("KIN_DECK_ALARM", "180")))
            try:
                if same:
                    res["deck_check"] = dict(prev["check"], reused="deck identical to the previous build")
                    deck_err = prev["check"]["err"]
                    raise StopIteration
                times = np.linspace(0, t_end, int(round(t_end / samp)) + 1)
                sim = deck_sim.simulate(run_dir, times)
                deck_err, dsp, _, _, _ = compare(rr, conv, {"Time(s)": times, **sim}, s_start, ev_times, samp)
                res["deck_check"] = {"err": deck_err, "worst_species": dsp}
            except StopIteration:
                pass
            except Exception as ex:
                res["deck_check"] = {"error": str(ex)[:200]}
            finally:
                signal.alarm(0)
        nuser = res["n_user_lines"]
        if att.get("complete") and att.get("n_compared", 1) == 0:
            status = "NOTHING_COMPARED"
        elif math.isfinite(score) and score < PASS:
            status = "PASS"
        elif math.isfinite(score) and score < APPROX:
            status = "APPROX"
        elif deck_err is not None and deck_err <= DECK_TOL:
            status = "DECK_VERIFIED_KINTECUS_NUMERICS"
        elif rr is None:
            status = "RAN_UNVERIFIED" if att.get("complete") else "KINTECUS_FAILED"
        elif not att.get("complete"):
            status = "KINTECUS_FAILED"
        else:
            status = "MISMATCH"
        if status == "DECK_VERIFIED_USER_BUG":
            conv.flag("partial", f"{nuser} USER rate lines: the deck reproduces the SBML (checked independently), but "
                                 "Kintecus V2027 evaluates only the last USER line correctly, so Kintecus output "
                                 "will differ until that is fixed")
        if blockers and status not in ("NOT_CONVERTIBLE", "ERROR"):
            status += "_PARTIAL"
        res["status"] = status
        res["switches"] = sw
        res["error"] = None if score == float("inf") else score
        # ---- kintpack
        desc = [f"{name}", f"BioModels {mid} (curated). SBML file: {MASTERS[mid]}.",
                f"Simulation: 0 to {K.num(t_end)} time units (SED-ML: {sedml or 'none'}"
                f"{'; duration chosen from the dynamics: ' + res['horizon_rule'] if res.get('horizon_rule') else ''}), "
                f"sampling {K.num(samp)}.",
                f"Run switches: {sw}",
                f"Verification against libRoadRunner (SBML reference): {status}" +
                (f", max |difference| = {score:.2e} of species peak (worst: {att.get('worst_species')})"
                 if math.isfinite(score) else ""),
                "SBML time units are used unchanged as Kintecus seconds."]
        if deck_err is not None:
            desc.insert(5, f"Independent check of this deck (SciPy integration with Kintecus semantics) against the "
                           f"SBML: max |difference| = {deck_err:.2e} of species peak.")
        if conv.flags:
            desc.append("Conversion notes:")
            desc += [f"  [{s}] {t}" for s, t in conv.flags]
        alias = [f"{r[0]} = {r[5]}" for r in conv.species_rows if r[0] != r[5] and not r[5].startswith("__")]
        if alias:
            desc.append("Renamed species (Kintecus name = SBML id): " + "; ".join(alias))
        pack = K.kintpack(model, species, parm, conv.ofiles, sw, "\n".join(desc),
                          datetime.datetime.now(datetime.timezone.utc).isoformat())
        os.makedirs(PACKS, exist_ok=True)
        json.dump(pack, open(f"{PACKS}/{mid}.kintpack", "w"), indent=1)
        res["kintpack"] = f"{mid}.kintpack"
    except K.Unsupported as ex:
        res["status"] = "NOT_CONVERTIBLE"
        res["flags"].append(["blocker", str(ex)])
    except Exception as ex:
        res["status"] = "ERROR"
        res["traceback"] = traceback.format_exc()[-1500:]
    res["seconds"] = round(time.time() - t0, 1)
    json.dump(res, open(f"{d}/result.json", "w"), indent=1, default=float)
    print(f"{mid} {res['status']} err={res.get('error')} {res['seconds']}s", flush=True)


K_EXE = "/home/claude/kin/Kintecus.exe"
if __name__ == "__main__":
    main(sys.argv[1])
