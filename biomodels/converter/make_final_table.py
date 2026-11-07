"""make_final_table.py -- FINAL_RESULTS_TABLE.csv / .txt: one row per curated BioModels entry with the Kintecus match
(% of species peak), SUCCESS (under 4%, James's rule of 2026-10-30; within 0.1% also shown), and a plain-language note on why it failed or how far off it is."""
import csv, textwrap, glob as _glob, glob, json, math, os, collections, datetime, textwrap

B = "/home/claude/kin/batch"
PASS_TOL = 4e-2        # SUCCESS = Kintecus under 4% of the SBML (James, 2026-10-30)
TOL = 1e-3             # strict match / deck verification (independent SciPy run of the pack)
CAP = 300
FIX_TS = os.path.getmtime(f"{B}/sbml2kin_batch_v2.py")   # converter backup made just before the #39/#999 fixes


def pct(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    if not math.isfinite(x):
        return "inf"
    p = x * 100
    return f"{p:.2g}" if p < 0.01 else f"{p:.3g}" if p < 100 else f"{p:.0f}"


def folder(r):
    st = r["status"].replace("_PARTIAL", "")
    if not r.get("kintpack"):
        return "(no pack)"
    if r["status"].startswith(("PASS", "NOTHING_COMPARED")):
        return "1_success_under_4pct"
    if st == "DECK_VERIFIED_KINTECUS_NUMERICS":
        return "2_deck_verified_kintecus_too_slow_or_differs"
    return "3_off_partial_or_unverified"


DEFTOL = sorted((os.path.basename(f)[5:-4].lstrip("0") for f in _glob.glob("/home/claude/kin/batch/references/*.csv")
                 if "default tolerances" in "".join(open(f).readlines()[:4])), key=int)


def short_int(sw):
    return "-INT:1" if "-INT:1" in sw else "-INT:5:10:0" if "-INT:5:10:0" in sw else "-INT:5" if "-INT:5" in sw else "-INT:2" if "-INT:2" in sw else sw


def attempt_text(a):
    sw = short_int(a.get("switches", ""))
    if str(a.get("where", "")).startswith("James's PC"):
        sw += " on James's PC"
    elif a.get("converter"):
        sw += " (deck before 11-01)"
    if a.get("fatal"):
        return f"{sw}: FATAL #{a['fatal']}"
    if a.get("nan"):
        return f"{sw}: reached the end time, but the output contains NaN"
    wd = a.get("watchdog")
    if wd:
        wd = wd.replace("stalled (no output row for 15 s", "stalled (no progress for 15 s")
        return f"{sw}: {wd}"
    if a.get("complete"):
        e = a.get("err")
        if e is not None and not math.isfinite(e):
            return f"{sw}: finished, but the output contains NaN/Inf"
        return f"{sw}: finished, off by {pct(e)}%" if e is not None else f"{sw}: finished"
    return f"{sw}: ended early at t={a.get('t_reached', 0):g}"


def why(r):
    st = r["status"]; base = st.replace("_PARTIAL", "")
    e = r.get("error"); dc = r.get("deck_check") or {}
    de = dc.get("err")
    best = ([r["attempts"][r["best_run"] - 1]] if r.get("best_run") and r["best_run"] <= len(r.get("attempts", [])) else
            [a for a in r.get("attempts", []) if a.get("switches") == r.get("switches")])
    best = best[0] if best else {}
    wsp = best.get("worst_species")
    att = "; ".join(attempt_text(a) for a in r.get("attempts", []))
    blockers = [t for s, t in r.get("flags", []) if s == "blocker"]
    partial = [t for s, t in r.get("flags", []) if s == "partial"]
    out = []
    if base == "PASS":
        if e is not None and e <= TOL:
            out.append(f"Kintecus matches the SBML within 0.1% (max {pct(e)}% of peak, worst species {wsp}).")
        else:
            out.append(f"PASS under the 4% rule: Kintecus is off by {pct(e)}% of peak (worst species {wsp}); "
                       f"before 2026-10-30 (0.1% rule) this was {r.get('status_strict_0.1pct', '?')}.")
            if de is not None and de <= TOL:
                out.append(f"The deck itself reproduces the SBML to {pct(de)}% (independent SciPy run of the same deck), "
                           f"so the remaining difference arises inside Kintecus (integration or evaluation).")
            elif de is not None:
                out.append(f"The deck itself is also off ({pct(de)}% in the independent run)" +
                           (", consistent with the unconverted parts below." if blockers else "."))
            elif dc.get("error"):
                out.append(f"Independent deck check not available ({dc['error'][:80]}).")
            if r["id"] in ("BIOMD0000000077", "BIOMD0000000281", "BIOMD0000000337", "BIOMD0000001019", "BIOMD0000001020"):
                out.append("Cause of the difference: Kintecus re-applies the last -SET value after the profile ends "
                           "(KINTECUS_ISSUES item 15).")
    elif base in ("APPROX", "MISMATCH"):
        out.append(f"Kintecus finishes but is off by {pct(e)}% of peak (worst species {wsp}).")
        if de is not None and de <= TOL:
            out.append(f"The deck itself reproduces the SBML to {pct(de)}% (independent SciPy run of the same deck), "
                       f"so the difference arises inside Kintecus (integration or evaluation).")
        elif de is not None:
            out.append(f"The deck itself is also off ({pct(de)}% in the independent run)" +
                       (", consistent with the unconverted parts below." if blockers else
                        ", so the conversion is the likely cause."))
        elif dc.get("error"):
            out.append(f"Independent deck check not available ({dc['error'][:80]}).")
    elif base == "DECK_VERIFIED_KINTECUS_NUMERICS" and any(a.get("complete") for a in r.get("attempts", [])):
        out.append(f"The deck reproduces the SBML to {pct(de)}% (independent SciPy run of the same deck), but "
                   f"Kintecus's output differs: {att}. The difference arises inside Kintecus (integration or "
                   f"evaluation).")
    elif base == "DECK_VERIFIED_KINTECUS_NUMERICS":
        out.append(f"The deck reproduces the SBML to {pct(de)}% (independent SciPy run of the same deck), but Kintecus "
                   f"did not finish within the time limit: {att}.")
    elif base == "KINTECUS_FAILED":
        out.append(f"Kintecus did not finish within the time limit ({att}).")
        if de is not None:
            out.append(f"Independent run of the deck: off by {pct(de)}%" +
                       (" (unconverted parts below)." if blockers else "."))
        elif dc.get("error"):
            out.append(f"The independent deck run also failed ({dc['error'][:80]}).")
    elif base == "NOTHING_COMPARED":
        out.append("Kintecus runs, but no displayed species changes over the run, so there is nothing to compare.")
    elif base == "RAN_UNVERIFIED":
        rr = [t for t in partial if "libRoadRunner" in t]
        out.append("Kintecus runs, but the SBML reference could not be simulated by libRoadRunner"
                   + (f" ({rr[0].split('(', 1)[-1].rstrip(')')[:90]})" if rr else "") + ", so it is unverified.")
    elif base == "NOT_CONVERTIBLE":
        out.append("No pack: " + ((r.get("error_text") or (blockers[0] if blockers else "")) or
                                  "nothing in the model maps onto Kintecus")[:160])
    elif base == "TIMEOUT_TOTAL":
        out.append("The conversion/verification pipeline exceeded its overall time limit for this model.")
    else:
        tb = (r.get("traceback") or "").strip().splitlines()
        out.append("Converter error: " + (tb[-1][:150] if tb else str(r.get("error"))[:150]))
    if blockers:
        out.append(f"Not converted ({len(blockers)}): " + "; ".join(b[:110] for b in blockers[:3]) +
                   (" ..." if len(blockers) > 3 else ""))
    mis = [t for t in partial if "mis-reads this long USER" in t]
    if mis:
        out.append(f"{len(mis)} long USER expression(s) are still mis-read by Kintecus V2027 (KINTECUS_ISSUES #2).")
    if r.get("horizon_rule"):
        out.append(f"Duration {r['t_end']:g} chosen from the dynamics (no curated SED-ML time course).")
    return " ".join(out)


def load():
    res = {}
    for f in glob.glob(f"{B}/out/*/result.json"):
        r = json.load(open(f)); res[r["id"]] = r
    return res


def build():
    res = load()
    rows = []
    for mid in sorted(res):
        r = res[mid]
        e = r.get("error"); de = (r.get("deck_check") or {}).get("err")
        blockers = [t for s, t in r.get("flags", []) if s == "blocker"]
        rows.append({
            "id": mid, "name": r.get("name", ""),
            "success_under_4pct": "YES" if r["status"].startswith(("PASS", "NOTHING_COMPARED")) else "NO",
            "within_0.1pct": "YES" if (e is not None and e <= TOL) else "NO",
            "kintecus_match_error_percent": pct(e),
            "deck_match_error_percent": pct(de),
            "status": r["status"], "folder": folder(r),
            "notes": why(r),
            "kintecus_runs": "; ".join(attempt_text(a) for a in r.get("attempts", [])),
            "not_converted_count": len(blockers),
            "duration": "" if r.get("t_end") is None else f"{r['t_end']:g}",
            "duration_source": "" if not r.get("kintpack") else ("dynamics" if r.get("horizon_rule") else "SED-ML"),
            "sbml_reactions": r.get("n_reactions_sbml", ""), "kintecus_lines": r.get("n_lines", ""),
            "species": r.get("n_species", ""),
            "built_with_latest_fixes": "yes" if os.path.getmtime(f"{B}/out/{mid}/result.json") >= FIX_TS and
                                        "previous build" not in " ".join(t for _, t in r.get("flags", [])) else "no",
            "kintpack": r.get("kintpack", ""),
        })
    with open(f"{B}/FINAL_RESULTS_TABLE.csv", "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    # text version
    c = collections.Counter(r["status"].replace("_PARTIAL", "") for r in res.values())
    ok = sum(r["success_under_4pct"] == "YES" for r in rows)
    ok01 = sum(r["within_0.1pct"] == "YES" for r in rows)
    packs = sum(1 for r in rows if r["kintpack"])
    L = [f"FINAL RESULTS - {len(rows)} curated BioModels entries -> {packs} .kintpack files   ({datetime.date.today()})",
         "", "SUCCESS = Kintecus output agrees with the original SBML (simulated by libRoadRunner) to under 4% (rule set",
         "          2026-10-30; it was 0.1% before - column within_0.1pct keeps the strict result):",
         "          max over the run of |Kintecus - SBML| / that species' peak, for every displayed species.",
         "Kintecus V2027: 20 packs whose events depend only on time (windows, periodic doses, clocks) were unrolled into",
         "time pairs and rerun on 2026-11-05; 152 packs with dosing/switch profiles were rerun on the 2026-11-05 08:40:46 build with",
         "Time/value-pair Constant Files (that build is needed for them).",
         "72 packs whose reaction lines grew with the 2,153-character limit were rerun on the 2026-11-03",
         "16:34:03 build (needed for the 51 packs with lines over 897 characters); the other packs with USER rate lines",
         "were last run on the 2026-10-30 22:51:25 build, the rest on",
         f"the 2026-10-26 15:30 build. Integrators tried in order: -ig:mass -INT:1 -obeymaxint, then "
         f"-INT:5 (DASPK), then for packs with signed",
         f"lines -INT:5:10:0 (DASPK without non-negativity), then for packs that still did not finish -INT:2 (Cash-Karp).",
         f"{CAP} s limit per attempt; -INT:2 got 600 s without the projected-time cut-off on the packs closest to finishing.",
         ] + textwrap.wrap(f"{len(DEFTOL)} references ({', '.join(DEFTOL)}) use libRoadRunner default tolerances because "
                           "CVODE fails at rtol 1e-10.", 118) + ["",
         f"  SUCCESS (under 4%)                                 {ok:5d}",
         f"    of which within 0.1%                             {ok01:5d}",
         f"  off by 4% or more                                  {c['MISMATCH']:5d}",
         f"  deck correct, Kintecus too slow / stops            {c['DECK_VERIFIED_KINTECUS_NUMERICS']:5d}",
         f"  Kintecus too slow / stops, deck not confirmed      {c['KINTECUS_FAILED']:5d}",
         f"  runs, SBML reference not simulable                 {c['RAN_UNVERIFIED']:5d}",
         f"  runs, nothing changes to compare                   {c['NOTHING_COMPARED']:5d}",
         f"  no pack (nothing convertible / converter error)    {c['NOT_CONVERTIBLE'] + c['ERROR'] + c['TIMEOUT_TOTAL']:5d}",
         f"  (of all packs, {sum(1 for r in res.values() if r['status'].endswith('_PARTIAL'))} have parts that could not "
         f"be converted - listed in the notes)",
         f"  {sum(1 for r in rows if r['built_with_latest_fixes'] == 'yes')} entries were rebuilt with the latest converter "
         f"(column built_with_latest_fixes); the rest keep their earlier build and result.", "",
         "Columns in FINAL_RESULTS_TABLE.csv: id, name, success_under_4pct, within_0.1pct, kintecus_match_error_percent (Kintecus vs",
         "SBML, % of peak), deck_match_error_percent (the same deck integrated independently in SciPy - shows whether the",
         "pack itself is right), status, folder, notes (why it failed / how far off), kintecus_runs (what each",
         "integrator attempt did), not_converted_count, duration, ...", "", "=" * 110]
    for r in rows:
        L.append(f"{r['id']}  {r['success_under_4pct']:3s}  Kintecus {(r['kintecus_match_error_percent'] or '-') + '%':>9s}"
                 f"  deck {(r['deck_match_error_percent'] or '-') + '%':>9s}  {r['status']}")
        L.append("    " + r["name"][:100])
        for ln in textwrap.wrap(r["notes"], 104):
            L.append("    " + ln)
    open(f"{B}/FINAL_RESULTS_TABLE.txt", "w", encoding="utf-8").write("\n".join(L) + "\n")
    return rows, c, ok


if __name__ == "__main__":
    rows, c, ok = build()
    print(open(f"{B}/FINAL_RESULTS_TABLE.txt").read()[:1500])
