Kintecus build 2026-11-03 12:39:17 (md5 e62fb158f2229e192a3725c36c074c13): checks of the two new features
Run each folder with:  Kintecus -ig:mass -INT:1 -obeymaxint -i

1. Time/value pairs: only the FIRST pair is used
   TIMEPAIRS_ADD_3PAIRS/     S has S.TXT-ADD; S.TXT = Time(s)/Concentrations, 100 1, 200 1, 300 1, END; 0 to 400 s.
                             Expected S = 1, 2, 3 after 100, 200, 300.  Seen: S = 1 after 100 and stays 1.
                             Same with -SET (100 -> 5, 300 -> 2: only the 100 s value is set), with times written
                             100.0, and without the END line. The first pair itself is right: applied exactly at its
                             time (output row at 100 still shows the old value, 101 the new one), correct amount.
   TIMEPAIRS_PINNED_2PAIRS/  P has P.TXT (no suffix) and initial value 7; P.TXT = 100 1, 200 3.
                             Expected P = 7 until 100, 1 until 200, then 3.  Seen: P = -987654000 from 0 to 100,
                             then 1 to the end. A pinned pair file whose first pair is at t = 0 (0 0, 100 1, 200 0)
                             stays at 0.
   A file with spaces instead of tabs ("Time(s)  Concentrations" / "100  1") gives FATAL Error #999.

2. Line length: the reaction field is still limited to 897 characters
   REACTION_FIELD_897/  source+S[USER;cA:(cA*(1+0.001+...))]==>X, reaction field 897 characters -> runs.
   REACTION_FIELD_898/  the same, 898 characters -> FATAL Error #24.
   A 2,100-character comment after a short reaction runs, so the whole row grew but the reaction-equation
   field did not. Same #24 for a 40-variable USER line of 1,387 characters (expression 689).

--- added after the 13:35:56 build (md5 baac5214596342089f774b79509785b1) ---
3. Reaction field: FIXED in 13:35:56. Up to 2153 characters (2159 with "1<TAB>0<TAB>0<TAB>") runs with correct
   values; 2154 -> FATAL #24.
4. USER rate evaluation is about 20x slower in both 11-03 builds (pure mass-action decks are not affected).
   SPEED_USER_BIOMD0000000010/ (10 USER lines, 0 to t_end, -INT:1):
       10-30 22:51 build: 2.1 s      11-03 12:39: 40.5 s      11-03 13:35: 40.0 s
   BIOMD0000000100 (8 USER lines): 2.0 s -> 39.6 s. BIOMD0000001028 (41 USER lines): t = 5.42 after 30 s on
   22:51, t = 0.002 after 30 s on 13:35. Looks like per-evaluation work that scales with the new 2,159-character
   buffers (copying, blank padding, TRIM/INDEX over the full length) rather than with the expression length.
5. A single pair is sometimes not applied at all:
   TIMEPAIRS_RUN_60S/  -ADD pair 14 s / 1.0, run 0-60 s, sampling 1 -> never applied.
   TIMEPAIRS_RUN_70S/  identical except a 70 s run -> applied at 14 s (seen on the 15 s output row).
   (60 s with sampling 0.1 or 30 s with sampling 0.01: applied.)
   TIMEPAIRS_BIOMD0000000439/  real pack: dose of L (0.1 at 14) as a pair, run 0-300 s, sampling 0.2857 ->
   never applied (L stays 1e-55). The same dose in the old spacing format is applied.
