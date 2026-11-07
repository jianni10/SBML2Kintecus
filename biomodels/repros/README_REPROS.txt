Kintecus V2027 bug repros. Rechecked with the build compiled 2026-10-30 22:51:25 (md5 1c592cfc...).

0c. (FIXED in 22:51 build) Variables per USER line: now 92 run; 93+ -> "*ERROR* Your User Defined function has
   improperly formatted line. You are missing a variable list" (a clean error, no memory corruption).
   USER_VARS_53_FATAL999 and USER_VARS_54_distinct_values now run. Checked: 52, 53, 54, 60, 90, 91, 92 run;
   93-98, 120, 200 give the error above. 60 variables in a 770-character list run, so the limit is the count.

0a. USER-line limit: now 2,043 in the 17:13 build                         (OK for all BioModels)
   2,043 USER lines run; 2,044 -> "SOFTWARE ERROR AT UFHS88 !". The largest BioModels pack has 1,961, so all fit.
   (943 in the 16:45 and 17:02 builds; 301 before that.)

0b. (FIXED in 22:51 build, see 0c) Variables per USER line: 52 in the 17:02 and 17:13 builds; 53+ corrupted memory
   USER_VARS_52_runs.kintpack            52 variables -> runs.
   USER_VARS_53_FATAL999.kintpack        53 variables -> "*ERROR9* There was a problem initializing the string
                                          function evaluator", then FATAL #999. In the variable list printed in
                                          VIEW.TXT, variable 1's name is 8 garbage bytes 9A 99 99 99 99 99 B9 3F
                                          = REAL*8 0.1.
   USER_VARS_54_distinct_values.kintpack  54 variables with distinct values 0.200..0.252: variable 1's name
                                          decodes to 0.251 = the value of variable 53.
   So the value of variable 53 is stored on top of the name of variable 1: a REAL*8 value array dimensioned 52
   (or a loop/limit of 52) sits just before the names in memory. BioModels affected: BIOMD0000000457 (90 variables).

Run each with:  Kintecus -ig:mass -INT:1 -obeymaxint

1. (FIXED in 16:45 build) USER-line limit: more than 301 USER lines in the 10-26 and 10-30 13:39 builds -> "SOFTWARE ERROR AT UFHS88 !"
   USER_LIMIT_301_runs.kintpack     301 USER lines (xa ==> xb, rate k_i*cxa)  -> runs to the end.
   USER_LIMIT_302_UFHS88.kintpack   the same plus one more USER line           -> console shows
                                    "SOFTWARE ERROR AT UFHS88 !" after "Total Reactions", nothing is
                                    integrated, KDONE.OUT says "FATAL ERROR ENCOUNTERED!!", no KERR.OUT.
   Only USER lines count: 301 USER lines + 200 mass-action lines runs; 200 mass-action + 302 USER fails.
   Identical or distinct USER expressions, 1 or 2 variables, product or NULL: the limit is always 301.
   Looks like a USER-function array dimensioned about 300.

   CORRECTION to my earlier "NULL bug" report: FATAL #999 "species name ... HAS NOT BEEN FOUND" at the 370th
   USER line is the same overflow showing up earlier, during parsing. My "control" with 368 filler lines did
   not run normally - it stopped with UFHS88 too (my script only checked that KDONE.OUT existed). With the
   3-variable last line: 300 fillers -> runs; 301 fillers (302 USER lines) -> UFHS88; 369 -> FATAL #999.
   NULL_BUG_repro_369.kintpack is kept as an example of the FATAL #999 form.

   BioModels packs that hit this (302 to 1,961 USER lines): BIOMD0000000235, 255, 332, 333, 334, 469, 470,
   471, 472, 473, 496, 497, 603, 604, 605, 606, 607, 608, 610. Raising the limit (to about 2,000) would let
   all 19 run.

2. USER expression buffer (~255 characters) with '^'                       (FIXED in the 10-30 build)
   CARET_BUG_fatal999_253chars.kintpack        one 253-character USER expression with 28 '^'
                                                10-26 build: FATAL #999.   10-30 build: prb(1 s) = 0.279129, correct.
   CARET_BUG_silently_wrong_289chars.kintpack  289 characters, 32 '^'
                                                10-26 build: prb(1 s) = 0.250514 (31% low, no warning).
                                                10-30 build: 0.361538, correct.
   CARET_BUG_control_same_rate_without_caret.kintpack  same rate written without '^' (417 characters)
                                                -> 0.361538 in both builds.

3. Constant File profile length                                             (2026-11-02, build 10-30 22:51)
   PROFILE_20002_VALUES/   A==>B (k = 1) for 10 s, species P pinned to PROF.TXT: spacing line -0.0005 and
                           20,002 values -> runs (CONC.TXT written).
   PROFILE_20003_VALUES/   the same with 20,003 values -> FATAL Error #94 while reading SPECIES.DAT.
   BioModels: BIOMD0000001029 needs 36,002 profile values (events on a fine time grid); 40,000 would cover it.
   The packs stop each profile just after its last event (negative spacing holds the last value).

4. -ADD profile entries applied one interval late far into a long profile          (2026-11-02, build 10-30 22:51)
   PROFILE_ADD_TIMING/     A==>B plus species D with ADD.TXT-ADD: spacing -0.1, 16,802 values, 1.0 at every
                           500th entry (t = 50, 100, ..., 1650 s); sampling 0.1, 0 to 1700 s, -INT:1.
   Expected: D steps up at the first output after k*0.1 (t = 50.1, 100.1, ...). Seen: 50.1, then 100.2 ... 300.2
   (one interval late), 350.1 ... 1200.1, then 1250.2 ... 1650.2 (late again). Looks like the profile clock is
   accumulated in limited precision. BIOMD0000000537 (sgpFc doses at 0.1, 336, ..., 1680 h, spacing 0.1): the
   1344 h dose lands at 1344.2 instead of 1344.1, which alone gives a 99% one-row difference; with those rows
   excluded the run is within 0.6% of libRoadRunner.
