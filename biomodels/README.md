# BioModels → Kintecus: 1,072 curated SBML models, converted and verified

Every curated [BioModels](https://www.ebi.ac.uk/biomodels/) entry converted to a Kintecus `.kintpack`, run in
Kintecus V2027 and compared with the original SBML simulated by [libRoadRunner](https://libroadrunner.org/).

| | Models |
|---|---|
| Curated entries | 1,072 |
| `.kintpack` files written | 1,045 |
| **SUCCESS** — Kintecus within 4% of the SBML | **942** |
| … of which within 0.1% | 921 |
| Deck reproduces the SBML (independent SciPy run), Kintecus too slow | 4 |
| Off, partial, or no SBML reference to check against | 99 |
| No pack (nothing convertible / converter error) | 27 |

Match error = max over species of |Kintecus − SBML| / peak of that species. Status as of 2026-11-06.

**Kintecus itself is not included.** It is copyrighted, proprietary software; get it from
[kintecus.com](http://www.kintecus.com/). The packs need the build of 2026-11-05 08:40:46 or later
(Time/value-pair Constant Files).

## Start here

* [`FINAL_RESULTS_TABLE.csv`](FINAL_RESULTS_TABLE.csv) — one row per model: status, match %, and why it failed
  or how far off it is. Plain-text version: [`FINAL_RESULTS_TABLE.txt`](FINAL_RESULTS_TABLE.txt).
* [`BioModels_Kintecus_Catalog.xlsx`](BioModels_Kintecus_Catalog.xlsx) — the same results with model names,
  descriptions and filters.
* [`plots/`](plots) — Kintecus vs. SBML overlay for every model that ran
  ([`plots/PLOT_INDEX.csv`](plots/PLOT_INDEX.csv)).
* [`README.txt`](README.txt) — the full release notes, conversion conventions and known Kintecus issues
  ([`KINTECUS_ISSUES.txt`](KINTECUS_ISSUES.txt)).

## Layout

| Folder / file | What it is |
|---|---|
| `1_success_under_4pct/` | packs that match the SBML to under 4% |
| `2_deck_verified_kintecus_too_slow_or_differs/` | deck verified, Kintecus too slow or stops |
| `3_off_partial_or_unverified/` | everything else that has a pack |
| `references/` | libRoadRunner reference time courses (CSV) used by `kintcheck.py` |
| `sbml2kintpack/` | **converter**: `sbml2kintpack.py` + `MAKE_KINTPACK.bat` (SBML → .kintpack); uses `sbml2kin_batch.py` and `release_settings.csv` |
| `converter/` | the batch pipeline that built and verified the collection (`run_model.py`, `gen_refs.py`, `deck_sim.py`, …) |
| `runkintpack/` | run a `.kintpack` in Kintecus from the command line |
| `kintpack2sbml/` | convert a `.kintpack` back to SBML |
| `kintcheck.py` | runs packs in Kintecus and compares with `references/` |
| `RUN_ALL_CHECKS.bat` | re-run every pack on your PC and compare with the SBML |
| `RUN_SLOW_JOBS.bat` | run the packs that were too slow here, without a time limit ([`SLOW_JOBS_README.txt`](SLOW_JOBS_README.txt)) |
| `RUN_REGRESSION.bat` | after a Kintecus change: re-run the passing packs and list any that now fail |
| `repros/` | small decks that reproduce Kintecus issues |

## Converting your own SBML

```bat
cd sbml2kintpack
MAKE_KINTPACK.bat path\to\model.xml
```

or `python sbml2kintpack.py model.xml --kintecus C:\Kintecus\Kintecus.exe`. Needs Python 3.8+ with
`python-libsbml`, `libroadrunner`, `sympy`, `numpy`. For a BioModels entry the duration, sampling, switches,
accuracy and deck variant of the verified pack are taken from `release_settings.csv`, so converting any of the
942 SUCCESS models from its original SBML gives the verified deck. Other SBML is converted the same way but is
not verified — check the pack description's conversion notes and "NOT converted" list.

## Checking a Kintecus build

```bat
RUN_REGRESSION.bat                         :: all passing packs
RUN_REGRESSION.bat --only BIOMD0000000535  :: selected models
```

Writes `regression_report_FAILED.csv`: packs that passed when built and now fail (FAILED) or pass with an error
more than 10× larger (WORSE).
