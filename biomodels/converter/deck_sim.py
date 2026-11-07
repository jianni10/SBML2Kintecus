"""
deck_sim.py -- integrate a generated Kintecus deck (MODEL/SPECIES/PARM + profile files) with SciPy, using the
documented/probed Kintecus semantics, and compare with libRoadRunner. Answers: is the DECK right, independent of
how Kintecus.exe evaluates it?

Semantics: line rate = (USER value or column-1 k) * prod(reactant concentrations); species change by
(product coeff - reactant coeff) * rate; a reactant coefficient n is also its order ([A]^n); column 7 numeric = pinned constant; profile file (first row = spacing,
negative = run once and hold): plain -> pinned value k on ((k-1)|s|, k|s|]; '-SET' -> value i set at t = i|s|;
'-ADD' -> value i added at t = i|s|.
"""
import json, math, os, re, sys
import numpy as np
from scipy.integrate import solve_ivp

FUN = {"exp": np.exp, "log": np.log, "log10": np.log10, "sqrt": np.sqrt, "abs": np.abs, "sin": np.sin,
       "cos": np.cos, "tan": np.tan, "asin": np.arcsin, "acos": np.arccos, "atn": np.arctan, "sinh": np.sinh,
       "cosh": np.cosh, "tanh": np.tanh, "floor": np.floor}


def rows(path):
    for ln in open(path, encoding="latin-1"):
        ln = ln.rstrip("\r\n")
        if not ln.strip() or ln.lstrip().startswith("#"):
            continue
        if ln.split("\t")[0].strip().upper() == "END":
            break
        yield ln.split("\t")


def side(s):
    out = []
    for part in s.split("+"):
        part = part.strip()
        if not part or part.upper() in ("NULL", "0 NULL"):
            continue
        m = re.match(r"^([0-9.]+(?:e-?\d+)?)?(.+)$", part)
        out.append((m.group(2), float(m.group(1)) if m.group(1) else 1.0))
    return out


def load_deck(d):
    names, init, disp, pins, profiles = [], {}, {}, {}, {}
    global OUTC
    OUTC = {}                                     # output-only species: name -> [(coefficient, species)]
    for r in rows(f"{d}/SPECIES.DAT"):
        nm = r[0].strip(); names.append(nm)
        init[nm] = float(r[2]); disp[nm] = r[3].strip().upper().startswith("Y")
        ssd = r[5].strip() if len(r) > 5 else ""
        m = re.match(r"(?i)^OUTC\[(.*)\]$", ssd)
        if m:
            terms = []
            for t in m.group(1).split(";"):
                t = t.strip(); sgn = -1.0 if t.startswith("-") else 1.0
                fac = [x for x in t.lstrip("-").split("*") if x]
                c, spn = sgn, None
                for x in fac:
                    try:
                        c *= float(x)
                    except ValueError:
                        spn = x
                terms.append((c, spn))
            OUTC[nm] = terms
        c7 = r[6].strip() if len(r) > 6 else "No"
        if c7.upper() in ("NO", "N", ""):
            continue
        try:
            pins[nm] = float(c7)
        except ValueError:
            mode = "SET" if c7.upper().endswith("-SET") else "ADD" if c7.upper().endswith("-ADD") else "PIN"
            fn = re.sub(r"-(SET|ADD)$", "", c7, flags=re.I)
            L = [l.rstrip("\r\n") for l in open(f"{d}/{fn}") if l.strip()]
            if "TIME" in L[0].upper():            # Time/value pairs (Kintecus build 2026-11-03), END-terminated
                pairs = []
                for l in L[1:]:
                    c = l.split("\t")
                    if c[0].strip().upper() == "END":
                        break
                    pairs.append((float(c[0]), float(c[1])))
                profiles[nm] = (mode, None, pairs)
            else:
                vals = [float(l.split("\t")[0]) for l in L]
                profiles[nm] = (mode, abs(vals[0]), vals[1:])
    lines = []
    for r in rows(f"{d}/MODEL.DAT"):
        k, rx = float(r[0]), r[3] if len(r) >= 5 else r[1]
        user = None
        m = re.search(r"\+?S\[USER;([^:]*):\((.*)\)\]", rx)
        if m:
            expr = m.group(2)
            rx = rx[:m.start()] + rx[m.end():]
            user = expr.replace("^", "**").lower()      # USER variables are case-insensitive
        lhs, rhs = rx.split("==>")
        lines.append((k, side(lhs), side(rhs), user))
    p = list(rows(f"{d}/PARM.DAT"))
    t_end = float(p[2][3])
    return names, init, disp, pins, profiles, lines, t_end


def simulate(d, times):
    names, init, disp, pins, profiles, lines, t_end = load_deck(d)
    idx = {n.upper(): i for i, n in enumerate(names)}
    y0 = np.array([init[n] for n in names], dtype=float)
    fixed = {}                                   # index -> current pinned value
    for n, v in pins.items():
        fixed[idx[n.upper()]] = v
    for n, (mode, s, vals) in profiles.items():
        if mode == "PIN":
            fixed[idx[n.upper()]] = vals[0] if s is not None else init[n]
    L = []
    for k, lhs, rhs, user in lines:
        net = {}
        for nm, c in lhs:
            net[idx[nm.upper()]] = net.get(idx[nm.upper()], 0) - c
        for nm, c in rhs:
            if nm.upper() in idx:
                net[idx[nm.upper()]] = net.get(idx[nm.upper()], 0) + c
        # reactant orders: integer coefficients as repeated factors, non-integer ones as powers (rate ~ [A]^n)
        L.append((k, [(idx[nm.upper()], c) for nm, c in lhs], net, user))
    # generate one Python function for the whole deck (much faster than eval per line)
    used = set()
    for _, _, _, user in L:
        if user is not None:
            used |= set(re.findall(r"\bc[a-z0-9_]+", user))
    var_of = {("c" + n).lower(): i for i, n in enumerate(names)}
    src = ["def _rhs(_t, y, dy):"]
    for v in sorted(used):
        if v in var_of:
            src.append(f"    {v} = y[{var_of[v]}]")
    src.append("    t = 298.0")
    for j, (k, reac, net, user) in enumerate(L):
        src.append(f"    r = {user if user is not None else repr(k)}")
        for i, c in reac:
            if c == int(c):
                src.extend([f"    r *= y[{i}]"] * int(c))
            else:
                src.append(f"    r *= y[{i}] ** {c!r}")
        for i, c in net.items():
            src.append(f"    dy[{i}] += {c!r} * r")
    src.append("    return dy")
    g = dict(FUN); g["__builtins__"] = _b
    exec(compile("\n".join(src), "<deck>", "exec"), g)
    _rhs = g["_rhs"]

    def rhs(t, y):
        for i, v in fixed.items():
            y[i] = v
        dy = np.zeros_like(y)
        with np.errstate(all="ignore"):
            _rhs(t, y, dy)
        for i in fixed:
            dy[i] = 0.0
        return dy

    # discontinuity schedule from profiles
    events = {}
    for n, (mode, s, vals) in profiles.items():
        i = idx[n.upper()]
        if s is None:                             # time/value pairs: each change once, at its time
            for t, v in vals:
                events.setdefault(round(t, 12), []).append((mode, i, v))
            continue
        if mode == "PIN":
            for kk in range(1, len(vals)):
                events.setdefault(round(kk * s, 12), []).append(("PIN", i, vals[kk]))
        else:
            for kk, v in enumerate(vals):
                events.setdefault(round(kk * s, 12), []).append((mode, i, v))
    y = y0.copy()
    for i, v in fixed.items():
        y[i] = v
    for mode, i, v in events.pop(0.0, []):
        if mode == "SET": y[i] = v
        elif mode == "ADD": y[i] += v
    bps = sorted(t for t in events if 0 < t < times[-1]) + [times[-1]]
    out = np.full((len(times), len(names)), np.nan)
    t0 = 0.0
    if times[0] == 0:
        out[0] = y
    for tb in bps:
        sel = (times > t0) & (times <= tb + 1e-12)
        if tb > t0:
            sol = solve_ivp(rhs, (t0, tb), y, method="LSODA", t_eval=np.clip(times[sel], t0, tb) if sel.any() else None,
                            rtol=1e-9, atol=1e-18)
            if sol.status < 0:
                raise RuntimeError(sol.message)
            if sel.any():
                out[sel] = sol.y.T
            y = sol.y[:, -1].copy() if sol.y.size else y
        for mode, i, v in events.get(tb, []):
            if mode == "PIN": fixed[i] = v
            elif mode == "SET": y[i] = v
            elif mode == "ADD": y[i] += v
        t0 = tb
    res = {n: out[:, j] for j, n in enumerate(names) if disp[n]}
    col = {n.upper(): j for j, n in enumerate(names)}
    for n, terms in OUTC.items():                 # OUTC[] species: output only (not part of the integration)
        if disp.get(n):
            res[n] = sum(c * (out[:, col[s_.upper()]] if s_ else 1.0) for c, s_ in terms)
    return res


OUTC = {}
import builtins as _b
_G = {"__builtins__": _b}


class _CI(dict):
    """case-insensitive variable lookup (Kintecus USER variables are case-insensitive)."""
    def __missing__(self, key):
        return dict.__getitem__(self, key.lower())


def check(mid):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import run_model as R
    res = json.load(open(f"{R.OUT}/{mid}/result.json"))
    att = ([res["attempts"][res["best_run"] - 1]] if res.get("best_run") and res["best_run"] <= len(res["attempts"]) else
           [a for a in res["attempts"] if a["switches"] == res.get("switches")])
    if not att:
        return None
    run_dir = f"{R.OUT}/{mid}/run_{res['attempts'].index(att[0]) + 1}"
    rr, init = R.rr_load(f"{R.ROOT}/{mid}/{R.MASTERS[mid]}")
    t_end = res["t_end"]; samp = res["sampling"]
    times = np.linspace(0, t_end, int(round(t_end / samp)) + 1)
    sim = simulate(run_dir, times)
    import sbml2kin_batch as K
    conv = K.convert(f"{R.ROOT}/{mid}/{R.MASTERS[mid]}", init, t_end)
    kin = {"Time(s)": times, **sim}
    ev = sorted({e["t"] for e in conv.events})
    err, wsp, per, n, nrows = R.compare(rr, conv, kin, res.get("t_start", 0.0), ev, samp)
    return err, wsp


if __name__ == "__main__":
    for mid in sys.argv[1:]:
        try:
            e = check(mid)
            print(mid, "deck-vs-SBML err", None if e is None else f"{e[0]:.3e} (worst {e[1]})", flush=True)
        except Exception as ex:
            print(mid, "deck sim failed:", str(ex)[:150], flush=True)
