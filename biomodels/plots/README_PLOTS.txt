KINTECUS vs SBML OVERLAY PLOTS - BioModels -> Kintecus kintpack conversion   (2026-10-31)
================================================================================
One folder per BioModels entry (biomodels_kintecus_vs_sbml_plots/<ID>/) with one plot:

  <ID>_overlay.png    Kintecus drawn on top of the original SBML:
                        SBML (libRoadRunner 2.11.0, rtol 1e-10, atol 1e-16; the same reference files
                        RUN_ALL_CHECKS.bat compares against) = solid colored lines
                        Kintecus (CONC.TXT of the run reported in FINAL_RESULTS_TABLE; switches in the title)
                        = black dashed lines over the matching colored line
                      The title gives the max difference as % of the species peak (the FINAL_RESULTS_TABLE score).

Species: up to 15 of the compared species (the displayed species of the pack), those with the largest peak
|concentration| in the SBML run - so up to 30 curves per plot (15 SBML + 15 Kintecus). When a model has more than
15, the title says "15 largest of N". One color per species (matplotlib tab10 plus gold, magenta, teal, rosy brown,
light green), in order of peak size; the legend lists them.
Time is in the model's time units (seconds for almost all BioModels).

No overlay plot is made when one of the two curves is missing (no Kintecus-only or SBML-only plots):
  - "no SBML reference"   libRoadRunner could not simulate the SBML (delays, fast reactions, qual models,
                          unreadable files), or nothing changes to compare - 22 entries
  - "no pack"             the model had nothing convertible - 28 entries
1,022 entries have an overlay plot.

PLOT_INDEX.csv lists every entry: success (under 4%), within 0.1% (earlier rule), Kintecus and deck (SciPy) match
error %, status, the overlay plot, the species plotted, and notes:
  - "Kintecus stopped at t=..."   the Kintecus run ended early (too slow / stalled); its dashed curves stop there
  - "not in CONC.TXT: ..."        a reference species missing from the Kintecus output
  - "no ... - no overlay plot"    see above

These plots are from the concentration runs of the conversion. They are unrelated to the -MECHV mechanism
validation, which writes no CONC.TXT.
