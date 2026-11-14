BIOMD0000000450 (converter of 2026-11-07): six switch species SW1-SW6, each a Time/value-pair file with 84-85 pairs
(daily windows), t_end 1000 s, sampling 1 s.

  Kintecus.exe -ig:mass -INT:1 -obeymaxint  -> CONC.TXT has rows at t = 0, 1, ..., 83 and then one row at 1001:
                                               nothing is written between t = 83 and the end (85 rows in all)
  Kintecus.exe -ig:mass -INT:5 -obeymaxint  -> a row every second (580 rows by t = 578 in 120 s here)

Same on the 2026-11-05 and 2026-11-12 builds. The values that are written match the SBML (0.0005%); comparing
by interpolating across the gap gives 0.59%. Probably the same -INT:1 / Time-pair handling as
PAIR_SKIP_INT1_BIOMD0000000554 (there the pair changes themselves are never applied).
