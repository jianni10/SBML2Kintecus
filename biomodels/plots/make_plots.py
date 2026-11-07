"""make_plots.py -- Kintecus over SBML (libRoadRunner) overlay plot for every BioModels pack.

For each model with both curves: plots/<ID>/<ID>_overlay.png (no Kintecus-only or SBML-only plots).
  Kintecus: CONC.TXT of the run reported in FINAL_RESULTS_TABLE (out/<ID>/run_N, N = the reported attempt).
  SBML:     references/<ID>.csv (libRoadRunner, rtol 1e-10, atol 1e-16), the same reference RUN_ALL_CHECKS uses.
  Species:  up to 15 compared species with the largest SBML peak |value|; each drawn twice (SBML solid colored,
            Kintecus dashed on top), so up to 30 curves per plot.
"""
import csv, json, os, sys, math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = "/home/claude/kin/batch"; OUT = "/home/claude/plots/biomodels_kintecus_vs_sbml_plots"
TOP = 15
# 15 distinguishable colors: matplotlib tab10 + 5 more (gold, magenta, teal, rosy brown, light green)
COLORS = list(plt.get_cmap("tab10").colors) + ["#e6c700", "#ff00ff", "#008080", "#bc8f8f", "#a6d854"]


def read_conc(path):
    with open(path, encoding="latin-1", errors="replace") as fh:
        lines = [l.rstrip("\r\n").split("\t") for l in fh if l.strip()]
    if not lines:
        return np.zeros(0), {}
    hdr = [h.strip() for h in lines[0]]
    idx = [(i, h) for i, h in enumerate(hdr) if h]
    rows = []
    for l in lines[1:]:
        try:
            v = [float(l[i]) for i, _ in idx]
        except (ValueError, IndexError):
            continue
        rows.append(v)
    a = np.array(rows) if rows else np.zeros((0, len(idx)))
    a[~np.isfinite(a)] = np.nan
    names = [h for _, h in idx]
    return a[:, 0], {n.upper(): (n, a[:, j]) for j, n in enumerate(names) if j > 0}


def read_ref(path):
    names, rows = None, []
    with open(path, newline="") as fh:
        for r in csv.reader(fh):
            if not r or r[0].startswith("#"):
                continue
            if names is None:
                names = r[1:]; continue
            rows.append([float(x) for x in r])
    a = np.array(rows)
    o = np.argsort(a[:, 0], kind="stable")
    a = a[o]
    return a[:, 0], {n.upper(): (n, a[:, j + 1]) for j, n in enumerate(names)}


def peak(y):
    y = y[np.isfinite(y)]
    return float(np.max(np.abs(y))) if y.size else 0.0


def style(ax, title, xlab):
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlab); ax.set_ylabel("Concentration")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)


def one(mid):
    try:
        return _one(mid)
    except Exception as e:
        plt.close("all")
        return {"id": mid, "name": "", "status": "", "overlay_plot": "",
                "species_plotted": "", "note": f"plot error: {type(e).__name__}: {e}"}


def _one(mid):
    r = json.load(open(f"{B}/out/{mid}/result.json"))
    name = r.get("name", "")[:90]
    info = {"id": mid, "name": r.get("name", ""), "status": r.get("status", ""), "overlay_plot": "",
            "species_plotted": "", "note": ""}
    if not r.get("kintpack"):
        info["note"] = "no pack (nothing convertible)"; return info
    a = r.get("attempts") or []
    sw = r.get("switches")
    run = ([r["best_run"] - 1] if r.get("best_run") and r["best_run"] <= len(a) else
           [i for i, x in enumerate(a) if x.get("switches") == sw])
    kin = None
    for i in (run + list(range(len(a) - 1, -1, -1))):
        p = f"{B}/out/{mid}/run_{i + 1}/CONC.TXT"
        if os.path.exists(p):
            tk, ck = read_conc(p)
            if tk.size:
                kin = (tk, ck, a[i].get("switches", "")); break
    ref = read_ref(f"{B}/references/{mid}.csv") if os.path.exists(f"{B}/references/{mid}.csv") else None
    if kin is None or ref is None:          # overlay only: both curves are needed
        info["note"] = "; ".join(x for x in (
            "" if kin else "no Kintecus output (run did not produce CONC.TXT)",
            "" if ref else "no SBML reference (libRoadRunner could not simulate it, or nothing to compare)") if x) + \
            " - no overlay plot"
        return info
    tk, ck, ksw = kin
    ts, cs = ref
    pool = {k: peak(v[1]) for k, v in cs.items()}          # the compared species, ranked by SBML peak |value|
    sel = sorted(pool, key=lambda k: -pool[k])[:TOP]
    disp = {k: cs[k][0] for k in sel}
    col = {k: COLORS[i % len(COLORS)] for i, k in enumerate(sel)}
    notes = []
    missing = [disp[k] for k in sel if k not in ck]
    if missing:
        notes.append("not in CONC.TXT: " + ", ".join(missing))
    t_end = r.get("t_end")
    stopped = t_end is not None and np.nanmax(tk) < t_end * (1 - 1e-6)
    if stopped:
        notes.append(f"Kintecus stopped at t={np.nanmax(tk):g} of {t_end:g}")
    err = r.get("error")
    mt = "" if err is None else f"   max diff {err * 100:.3g}% of peak" if math.isfinite(err) else "   (non-finite)"
    more = f" ({TOP} largest of {len(pool)})" if len(pool) > TOP else ""
    d = f"{OUT}/{mid}"; os.makedirs(d, exist_ok=True)
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(figsize=(10, 5.5))
    dark_bg = False
    for k in sel:
        ax.plot(ts, cs[k][1], color=col[k], lw=3.2, ls="-", label=disp[k])
        if k in ck:
            c3 = matplotlib.colors.to_rgb(col[k])
            lum = 0.2126 * c3[0] + 0.7152 * c3[1] + 0.0722 * c3[2]
            kc = "white" if lum < 0.15 else "black"        # white dashes on a very dark SBML line
            dark_bg |= kc == "white"
            ax.plot(tk, ck[k][1], color=kc, lw=1.2, ls=(0, (5, 3)))
    stop_t = f"\n[Kintecus stopped at t={np.nanmax(tk):g} of {t_end:g}]" if stopped else ""
    sub = "; ".join(x for x in (mt.strip(), more.strip(" ()")) if x)
    style(ax, f"{mid} - {name}\nSBML / libRoadRunner (solid, colored) vs Kintecus {ksw} (dashed)"
              + (f"\n{sub}" if sub else "") + stop_t, "Time (s)")
    h, l = ax.get_legend_handles_labels()
    h.append(Line2D([], [], color="black", lw=1.2, ls=(0, (5, 3)))); l.append("Kintecus (each species)")
    if dark_bg:
        h.append(Line2D([], [], color="white", lw=1.2, ls=(0, (5, 3)))); l.append("Kintecus on dark lines")
    ax.legend(h, l, fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    fig.tight_layout(); fig.savefig(f"{d}/{mid}_overlay.png", dpi=100); plt.close(fig)
    info["overlay_plot"] = f"{mid}/{mid}_overlay.png"
    info["species_plotted"] = "; ".join(disp[k] for k in sel)
    info["note"] = "; ".join(notes)
    return info


if __name__ == "__main__":
    ids = sys.argv[1:] or sorted(os.path.basename(os.path.dirname(f)) for f in
                                 __import__("glob").glob(f"{B}/out/*/result.json"))
    os.makedirs(OUT, exist_ok=True)
    res = []
    with ProcessPoolExecutor(2) as ex:
        for i, x in enumerate(ex.map(one, ids, chunksize=4)):
            res.append(x)
            if i % 100 == 0:
                print(i, x["id"], flush=True)
    if len(sys.argv) == 1:
        tab = {row["id"]: row for row in csv.DictReader(open(f"{B}/FINAL_RESULTS_TABLE.csv", encoding="utf-8-sig"))}
        for x in res:
            t = tab.get(x["id"], {})
            x["success_under_4pct"] = t.get("success_under_4pct", "")
            x["within_0.1pct"] = t.get("within_0.1pct", "")
            x["kintecus_match_error_percent"] = t.get("kintecus_match_error_percent", "")
            x["deck_match_error_percent"] = t.get("deck_match_error_percent", "")
            x["status"] = t.get("status", x.get("status", ""))
        f = ["id", "name", "success_under_4pct", "within_0.1pct", "kintecus_match_error_percent", "deck_match_error_percent", "status", "overlay_plot", "species_plotted", "note"]
        with open(f"{OUT}/PLOT_INDEX.csv", "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=f, extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("done", len(res))
