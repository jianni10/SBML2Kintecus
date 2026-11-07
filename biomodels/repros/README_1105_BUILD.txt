Kintecus build 2026-11-05 08:40:46 (md5 9b45994dfd259bf24ec1fffc648345ce): time/value-pair checks
All runs: Kintecus -ig:mass -INT:1 -obeymaxint -i; the same results with -INT:5 and -INT:2.

FIXED (all correct, exact values, applied at the pair time):
  -ADD, 3 pairs (100/200/300 s, 1 each)                       S = 1, 2, 3
  -SET, 2 pairs (100 -> 5, 300 -> 2), initial 0.5             0.5, 5, 2
  single -ADD pair at 14 s: runs 0-60 (sampling 1), 0-70, 0-300 (sampling 0.2857142858), 0-60 (0.1): applied
  -ADD pair at t = 0 and at t = 1e-12 (below the first step): applied at the start
  1,000 -ADD pairs (every 1 s): every output row right, final 1000
  1,000 pinned pairs: every output row right
  50 pair files in one model (two pairs each): all applied
  no END line, decimal times (100.0, 200.5): right
  reacting species (S ==> Y, k = 0.01) with two doses and one pair after t_end: within 5e-6 of the exact solution
  pinned pairs 0/100/200 (0, 1, 0): right
  TIMEPAIRS_ADD_3PAIRS, TIMEPAIRS_RUN_60S, TIMEPAIRS_BIOMD0000000439: now right

STILL OPEN (minor; the converter never writes pinned pair files):
  TIMEPAIRS_PINNED_2PAIRS: pinned P, initial value 7, pairs 100 -> 1, 200 -> 3. P = -987654000 until 100 s, then
  1, then 3. Expected 7 until 100 s (the SPECIES.DAT initial value), as -SET files do.

Not a pair issue (harness mistakes found while testing, noted so nobody repeats them):
  * a blank Const column in SPECIES.DAT next to a pinned species -> FATAL #11 (use "No").
  * a model whose only reaction consumes the pinned species ("P ==> Y", one dynamic species) ->
    "Software error, SPEC.LE.2"; same on the 11-03 16:34 build.

Regression: 29 passing packs (USER rates, long lines, spacing profiles) give the same match errors as the 11-03
16:34 build (to 4 digits); reaction field limit still 2,153.

New on 11-05 (question): the output grid restarts at every time pair.
  TIMEPAIRS_OUTPUT_GRID/pairs: A ==> B, one -ADD pair at t = 14, 60 s run (output step 0.9836). CONC.TXT rows:
  ... 12.7869, 13.7705, 14.0, 14.9836, 15.9672 ... 59.2459, 60.0: a row at 14.0, then the grid starts again
  from 14. TIMEPAIRS_OUTPUT_GRID/spacing (same dose as a spacing profile) keeps the grid: 13.7705, 14.7541, ...
  The values are correct either way; only the row times move. In the BioModels packs every row after the
  first pair time is off the PARM grid (BIOMD0000000535: 660.1, 672, 702, 732 instead of 660, 690, 720), so a
  comparison on the PARM grid interpolates across the dose (44% on one row). The kintcheck references for the
  pair packs are now written at the times Kintecus actually prints. If the grid restart is not intended,
  keeping the original grid and adding the row at the pair time would be the expected output.
