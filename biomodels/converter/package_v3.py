"""package_v3.py -- final deliverable zip: packs by result, FINAL_RESULTS_TABLE (csv+txt), README, DESCRIPTIONS,
KINTECUS_ISSUES, references/ + kintcheck.py + RUN_ALL_CHECKS.bat (self-check on the user's PC), converter source."""
import csv, glob, io, json, os, sys, zipfile, datetime, collections
sys.path.insert(0, "/home/claude/kin/batch")
import make_final_table as FT

B = "/home/claude/kin/batch"
rows, counts, n_ok = FT.build()
res = FT.load()
today = datetime.date.today()
cnt = collections.Counter(r["folder"] for r in rows)
packs = [r for r in rows if r["kintpack"]]
refs = {os.path.basename(f)[:-4] for f in glob.glob(f"{B}/references/*.csv")}
n_ref = sum(1 for r in packs if r["id"] in refs)
deck_ok = sum(1 for r in rows if r["kintpack"] and (r["success_under_4pct"] == "YES" or
                                                    r["status"].startswith("DECK_VERIFIED")))

README = f"""BioModels -> Kintecus .kintpack collection
==========================================
Built {today}. Source: the curated BioModels set as mirrored in github.com/sys-bio/temp-biomodels
(entries up to BIOMD0000001080, mirror dated June 2025).

{len(rows)} curated models -> {len(packs)} .kintpack files.
SUCCESS (Kintecus output under 4% of the original SBML; rule set 2026-10-30): {n_ok} packs
  ({sum(1 for r in rows if r["within_0.1pct"] == "YES")} of them within 0.1%, the earlier rule; column within_0.1pct in the table).
Packs shown to reproduce the SBML (SUCCESS, or folder 2 where an independent SciPy run of the deck matches
  within 0.1% but Kintecus is too slow or differs): {deck_ok}.
{sum(1 for r in rows if r["built_with_latest_fixes"] == "yes")} entries (including every former #39/#999 pack) were rebuilt with the latest converter
fixes; the others keep their earlier build and result (column built_with_latest_fixes in the table).

START HERE: FINAL_RESULTS_TABLE.csv (open in Excel) or FINAL_RESULTS_TABLE.txt - one row per model with the
match %, SUCCESS yes/no, and a note on why it failed or how far off it is.

WHAT IS IN THIS ZIP   (full table of contents: CONTENTS.txt)
  1_success_under_4pct/                          {cnt['1_success_under_4pct']:4d}  Kintecus output matches the SBML to under 4%
  2_deck_verified_kintecus_too_slow_or_differs/  {cnt['2_deck_verified_kintecus_too_slow_or_differs']:4d}  the deck reproduces the SBML (independent
                                                       check), but Kintecus is too slow, stops, or differs
  3_off_partial_or_unverified/                   {cnt['3_off_partial_or_unverified']:4d}  off by 4% or more, too slow, parts not convertible,
                                                       or no SBML reference to check against
  (no pack)                                      {cnt['(no pack)']:4d}  nothing convertible / converter error
  FINAL_RESULTS_TABLE.csv / .txt   the results table (see above)
  RUN_ALL_CHECKS.bat               re-run every pack on YOUR PC and compare with the SBML reference (below)
  RUN_SLOW_JOBS.bat                run the packs that did not finish here, with no time limit (SLOW_JOBS_README.txt)
  RUN_REGRESSION.bat               after a Kintecus change: re-run the passing packs and list the ones that now fail
                                   (regression_report_FAILED.csv; FAILED = over 4%, WORSE = error grew > 10x)
  kintcheck.py                     the checker the .bat runs (plain Python 3, no packages needed)
  (original SBML + SED-ML inputs: separate download biomodels_original_sbml.zip, one folder per entry)
  references/                      SBML reference output for {n_ref} packs (libRoadRunner), used by kintcheck.py
  DESCRIPTIONS.txt                 one entry per model: results, switches, what was not converted
  KINTECUS_ISSUES.txt              Kintecus V2027 issues found, with minimal repros
  results_detail.csv               the full per-model data behind the table
  sbml2kintpack/                   convert ANY SBML file to a .kintpack on your PC (see below)
  runkintpack/                     unpack any .kintpack and run Kintecus on it, with field overrides (see below)
  kintpack2sbml/                   export a .kintpack or Kintecus deck to SBML / SED-ML / Kintecus XML
  converter/                       the scripts that built and verified every pack

CHECKING THE PACKS ON YOUR OWN PC (RUN_ALL_CHECKS.bat)
  1. Unzip everything into one folder (keep the folder structure).
  2. Install Python 3.8+ from python.org if you don't have it (no extra packages are needed).
  3. Open RUN_ALL_CHECKS.bat in Notepad and set the Kintecus.exe path on the line
        if "%KINTECUS_EXE%"=="" set "KINTECUS_EXE=C:\\Kintecus\\Kintecus.exe"
     (or put Kintecus.exe next to the .bat, or answer the prompt when it asks).
  4. Double-click RUN_ALL_CHECKS.bat. For each pack it writes the DAT/profile files into
     kintcheck_work\\<ID>\\ exactly as the GUI does, runs Kintecus with the pack's switches plus -i, watches
     KDONE.OUT / KERR.OUT / CONC.TXT, and compares CONC.TXT with references\\<ID>.csv.
  5. Results: kintcheck_report.txt (summary) and kintcheck_report.csv (per pack: PASS/OFF/STALLED/...,
     match error %, worst species, Kintecus message, and the result this pack gave when it was built, so a
     new Kintecus version's changes stand out). Failed runs stay in kintcheck_work\\ for inspection.
  Options (add after the .bat name in a command prompt, or edit the .bat):
     --folder 1_success          only packs in matching folders     --only BIOMD0000000010 ...  only these IDs
     --stall 120 --cap 1800      more patience (seconds)            --tol 0.04                   threshold (4%)
     --switches "-ig:mass -INT:5:59:1e-20 -obeymaxint"  force other switches     --show  show Kintecus windows
  A full run of all packs takes hours (each stuck pack waits --stall seconds; slow ones up to --cap).
  Start with  RUN_ALL_CHECKS.bat --folder 1_success  (about {cnt['1_success_under_4pct']} packs).

FINISHING THE PACKS THAT RAN OUT OF TIME HERE (RUN_SLOW_JOBS.bat)
  SLOW_JOBS.csv lists every pack that ran out of time on all integrators tried here, with the integrator that
  progressed fastest and a run-time estimate. RUN_SLOW_JOBS.bat runs them with no time limit, fastest first, prints
  progress every minute and compares each finished run with its SBML reference (slow_jobs_report.csv). Full
  instructions: SLOW_JOBS_README.txt. Packs that failed on a Kintecus error instead are in NOT_TIME_PROBLEMS.csv.
  Note: kintcheck.py was tested here under Wine on Linux; the native-Windows launch (hidden console via
  CREATE_NEW_CONSOLE) could not be exercised here - if Kintecus windows misbehave, try --show.

CONVERTING YOUR OWN SBML FILES (sbml2kintpack/)
  Install once:  pip install python-libsbml sympy libroadrunner numpy
                 (libroadrunner + numpy are optional but recommended: exact initial values and automatic duration)
  Drag SBML files onto sbml2kintpack\\MAKE_KINTPACK.bat, or from a command prompt:
     python sbml2kintpack.py model.xml                         duration from a .sedml next to the model
     python sbml2kintpack.py model.xml --t-end 500 --points 1000
     python sbml2kintpack.py model.xml --kintecus "C:\\Kintecus\\Kintecus.exe"
  Output next to the model: <name>.kintpack (load in the GUI) and <name>_deck\\ (MODEL.DAT, SPECIES.DAT,
  PARM.DAT, profile files). --kintecus runs the long-USER-expression check in Kintecus itself and repairs
  lines Kintecus mis-reads (recommended; MAKE_KINTPACK.bat does it when it finds Kintecus.exe - set the
  path at the top of the .bat). For a BioModels entry (ID in the file or folder name, or the model's own
  annotation) the duration, sampling, switches, accuracy and deck variant (shifted states etc.) of the released
  pack are taken from sbml2kintpack\release_settings.csv: converting any of the 890 SUCCESS models from
  biomodels_original_sbml.zip this way was checked here to give the verified deck (838 identical; 52 written by
  the current converter instead of an older one, all 52 run in Kintecus and PASS). --no-release-settings skips
  the table. States that go below zero are written shifted with an OUTC[] output (--no-shift: signed lines).
  Other SBML files are converted the same way but are not verified - check the description's conversion notes
  and "NOT converted" list.

RUNNING A PACK FROM THE COMMAND LINE (runkintpack/)
  runkintpack.bat PACK.kintpack [options]        (or: python runkintpack.py ...; Python 3.8+, no packages)
  Writes every sheet and file of the pack - MODEL.DAT, SPECIES.DAT, PARM.DAT, THERM.DAT, TRANSPORT.DAT,
  fitdata and perturbation/profile files - tab-delimited into <pack>_run\\ (or -o DIR) and runs Kintecus.exe
  there. The .kintpack itself is never changed.
     -k "SWITCHES"   override the pack's Kintecus switches       -p PATH   Kintecus.exe (overrides the pack's)
     -i              international inputs: 4,34 -> 4.34 (--set fields are then separated by ';')
     --set species.H2O.initialconcentration=4.34      --set parm.temperature=897
     --set model.12="3.4e10,0,0,A+B==>C,my comment"  --set model.12.A=3.4e10
     --set fitdata001.5="0,1.25,3.5"                  --set therm.3.2=NASA     --set transport.4="..."
     --list (show settable fields)  --check (validate --set only)  -u (unpack only)
     --batch [--stall S] [--cap S]  (unattended: hidden console, stops stuck runs)
  runkintpack.bat --help lists everything; runkintpack\\RUNKINTPACK_MANUAL.txt has the full manual with examples.

EXPORTING TO SBML (kintpack2sbml/)
  kintpack2sbml.bat INPUT [options]   INPUT = a .kintpack, a Kintecus deck folder, or its MODEL.DAT
  Writes NAME.xml (SBML L3V2), NAME.sedml (time course) and NAME_kintecus_unsupported.txt (everything SBML
  cannot hold, with the original text; also in an XML comment and the model notes). --format both adds a
  lossless NAME.kintecus.xml. --validate runs the libSBML check. Needs: pip install python-libsbml.
  Full manual: kintpack2sbml\\KINTPACK2SBML_MANUAL.txt

HOW THE MATCH IS SCORED
  Reference = the original SBML simulated by libRoadRunner (relative tolerance 1e-10, absolute 1e-16) on the
  pack's own sampling grid. For every displayed species: max over the run of |Kintecus - SBML| divided by that
  species' peak |value|; the pack's match error is the worst species. SUCCESS < 4% (0.1% before 2026-10-30).
  Rows within half a sampling interval of an event time are skipped (both sides jump there), and so is the part
  before the SED-ML output start time. The deck column in the table is the same score for the deck integrated
  independently in SciPy with Kintecus semantics: if the deck matches but Kintecus does not, the difference
  arises inside Kintecus. The shipped references keep up to 1000 rows per pack (plus each species' extreme
  points); the table was scored on every row.

HOW TO USE A PACK (Kintecus GUI)
  Home page -> "Load Package..." (or drag the .kintpack onto the drop zone). All grids load, plus any O_ files
  (Constant File profiles for events/switches). Run switches are stored in the pack, normally
     -ig:mass -INT:1 -obeymaxint      or, where that failed,   -ig:mass -INT:5 -obeymaxint (DASPK)
     or, where both ran out of time,  -ig:mass -INT:2 -obeymaxint (Cash-Karp Runge-Kutta, no Jacobian)
  Read the pack's Description box first: verification results, time span, conversion notes, what was NOT
  converted, renamed species.

CONVERSION CONVENTIONS (also in each pack's description)
  * SBML time units are used unchanged as Kintecus seconds; concentrations as given (MOLES/LITER label).
  * Mass-action laws -> plain rate constants; other laws -> A+S[USER;cX,...:(expr)]==>B lines.
  * Each reaction's comment is on its own '#' line just above it. Reaction fields are at most 2,153 characters
    (longest in the collection: 2,148, BIOMD0000000406); 51 packs have fields over 897 and need Kintecus V2027
    compiled 2026-11-03 16:34:03 or later (older builds: 897, FATAL #24 above it).
  * At most 3 reactant entries per line (Kintecus limit, FATAL #39): rate modifiers are kept inside the USER
    expression when a line would exceed it, and reactions with 4+ consumed reactants are split over several
    lines carrying the same rate.
  * Long USER expressions (Kintecus mis-reads them, see KINTECUS_ISSUES #2): constants shortened to 7 digits,
    powers written without '^' (x*x, exp(y*log(x))), species in them given short names (Q1, Q2, ...; listed in
    the description and the SPECIES comments column), and each long line's value checked in Kintecus itself
    against the intended value, re-padded with redundant parentheses where Kintecus mis-read it.
  * Kintecus species names are at most 26 characters (a USER variable c<name> may have 27). Packs rebuilt
    for this release need Kintecus V2027 compiled 2026-10-30 22:51:25 or later: older builds limit USER
    variables to 10 characters and misread function names inside them (before 10-26 15:30; KINTECUS_ISSUES 4
    and 5), allow 301 or 943 USER lines (before 10-30 17:13; now 2,043) and 52 variables per USER line
    (before 10-30 22:51; now 92, and 93+ gives a clean "missing a variable list" error).
    Packs with reaction fields over 897 characters need the 2026-11-03 16:34:03 build (md5 59fca291...),
    which also carries the faster USER interpreters (KINTECUS_ISSUES 22). Packs whose Constant Files are
    Time/value pairs (first line contains "Time") need the 2026-11-05 08:40:46 build (md5 9b45994d...;
    KINTECUS_ISSUES 21). Kintecus holds 95 pair files of up to 1,052 pairs (TSX1/TSY1(95,1052)); the
    largest model uses 64 files (BIOMD0000000235) and the largest file 7 pairs.
  * Long rates (since 2026-11-06): a rate law the SBML writes in factored form ((1+a/Ka)(1+b/Kb)...) is kept in
    that form (sympy's cancel/expand used to multiply such laws out to 5,000-20,000 characters); a rate that is
    still longer than the 2,153-character field and is a sum (or a factor times a sum) is written as several lines
    with the same reactants and products, each carrying part of the rate (rates add, so this is exact).
  * A reaction id used inside another expression is that reaction's rate (its kinetic law); constants beyond
    the double range (e.g. exp(3800) inside a switch) are written as +-1e308, which is what IEEE arithmetic in
    Kintecus gives; a time trigger whose threshold the events themselves change (time >= t_int*cups + t_lag,
    cups = cups + 1) is unrolled like the time-only events below.
  * Step functions, signs and exp() (since 2026-11-07; why the slow jobs were slow). Kintecus USER has no
    conditionals, so SBML piecewise/relations on the state are written with H(z) = (z/(abs(z)+1e-300)+1)/2.
    That step, and max(z,0) = (z+abs(z))/2, are now kept as one factor while rates are simplified; before,
    sympy multiplied them out into terms of size ~1/|z| on separate lines that cancel, so near the switch
    point the rate was rounding noise and the integrators crawled (deck round trip through CVODE: 158,707 steps
    for BIOMD0000000603 against 1,767 for the SBML; now 1,832). Terms of opposite sign that depend on the same
    species stay one rate instead of separate (signed) lines; exp() is kept whole (sympy wrote exp(3800-4t) as
    1e308*exp(-4t), which overflows to inf/NaN); conditions on time alone (daily light windows with
    floor(time/24), time - c >= 0, i == k) become switch species with Time/value pairs; and a negative
    product of a sum is printed -1*(a+b), not -1*a+b. Packs kept from before this change are rebuilt the old
    way by sbml2kintpack (release_settings.csv variant "...+pre1107").
  * Events that depend only on time (since 2026-11-05): triggers built from time, constants and variables that
    only the events themselves change - time windows ((time > 50) && (time <= 300)), periodic doses written with
    % or floor ((time - t1) % P == 0), self-resetting clocks (Day = Day + 24), exact instants (time == 20) - are
    unrolled by an exact discrete simulation (piecewise-affine triggers: breakpoints of floor/% and the root of
    each relation inside each segment) into their firing times up to the end time, and written as -SET/-ADD
    Time/value pairs. 20 models; firing times checked against libRoadRunner. Events whose triggers use a
    continuously changing variable (thresholds on concentrations), delays, or values taken from the state at
    the firing time are still not converted (they need event handling inside Kintecus).
  * Non-integer reactant stoichiometry (since 2026-10-31): a coefficient up to 4 is written as a Kintecus
    coefficient ('0.25CaER+...': Kintecus uses [CaER]^0.25 in the rate and consumes 0.25 x rate), and the USER
    rate is divided by [CaER]^0.25 to match. A larger one (e.g. 2923.98) is consumed on its own line
    ('X+S[USER;...:(n*rate/cX)]==>sink') so [X]^n never over- or underflows.
  * Rate order different from the stoichiometry (since 2026-11-01). Kintecus ties a reactant's order to its
    coefficient; when the SBML rate is not proportional to [A]^n, A is consumed on its own line:
    first order in A -> 'A+S[USER;n*rate/cA]==>sink'; zero order in A -> 'source+S[USER;-n*rate]==>A' (a signed
    line, exactly as the SBML, so A may go below zero where the SBML's does; such packs also get -INT:5:10:0).
    Rates where no term is proportional (lin-log kinetics) are written once, signed. This replaced the
    divided form rate*[A]/([A]+1e-30), whose on/off switch at [A] ~ 1e-30 stalled the integrators, and the
    huge powers ((cA+1e-30)^99) behind FATAL #167. Time conditions inside SBML assignment rules now become
    0/1 switch species (Constant Files) like the others, not abs() step functions.
    Where the new deck does worse, run_model also tries the deck written before 11-01 and keeps the better one.
  * Variables the SBML takes below zero (membrane potentials, deviations from basal, species of models whose
    published rates let them go negative; detected from the libRoadRunner run, more than 0.1% of the peak):
    Kintecus integrates X_o = X + shift (shift = 1-2-5 step >= 2 x |min|, so X_o stays positive), rates use
    (cX_o-shift), and the output-only species X = OUTC[X_o;-shift*source] gives the value itself (CONC.TXT, the
    plots and the references use X). The SPECIES comments give the shift. If X is a concentration, the
    published model is defective (it lets X go negative).
  * No NULL: an inflow is written source==>X and an outflow X==>sink, with source and sink held constant at
    1.0 in SPECIES.DAT (an SBML species named Source or Sink becomes Source1 / Sink1).
  * PARM.DAT sampling interval and maximum integration time are rounded up in the 10th digit, so Constant
    File changes land exactly on their event times (KINTECUS_ISSUES 12).
  * Species never changed by a reaction are pinned with SPECIES.DAT column 7.
  * Time events -> Constant File profiles (pinned, -SET or -ADD); conditions on time -> 0/1 switch species.
    State-triggered events, delays, algebraic rules and variable compartments are NOT converted and are
    listed in the notes.
  * Entries whose SED-ML was only an auto-generated placeholder got a duration chosen from the dynamics.
"""

ISSUES = """Kintecus V2027 issues found while converting BioModels (minimal repros; all deterministic)
=================================================================================

FIXED in the build compiled 2026-10-25 12:23
1. Species variables in a USER line read the wrong concentrations unless that line is the LAST USER line.
   1  B+S[USER;cA:(0.1*cA)]==>NULL         A=2, D=1 pinned (SPECIES col 7); B=C=1; 10 s
   1  C+S[USER;cD:(0.3*cD)]==>NULL
   Original build: B(10)=0.7285. Fixed build: B(10)=0.135335 (exact).

FIXED in the build compiled 2026-10-26 14:44 (and kept in 15:30)
5. Function and constant names were recognized INSIDE USER variable names: 'pi' was substituted
   (cnfatpinu1 -> cnfat3.14159265nu1), and a variable containing ln, exp, log, sqrt, abs, sin, cos, tan or
   floor gave FATAL #999 (A+S[USER;cCln:(0.2*cCln)]==>Z). Now ln, exp, log, log10, sqrt, abs, sin, cos, tan,
   floor, asin, acos, sinh, atn and pi all work at the start, middle or end of a name, also followed by
   digits (cPI3K, cLN2). The converter no longer renames such species (Glucose was Glucxose1 before).
   USER evaluation is also 3-7x faster (BIOMD0000000236: 154 s -> 23 s, BIOMD0000000015: 115 s -> 28 s);
   96 packs that were too slow before now finish and match the SBML.

FIXED in the build compiled 2026-10-26 15:30
4. USER variable names were limited to 10 characters ('c' + a 9-character species name; longer gave
   FATAL #999). Now variables up to 27 characters work - 'c' + a 26-character species name:
       1 0 0 A+S[USER;cXXXXXXXXXXXXXXXXXXXXXXXXX1:(0.2*cXXXXXXXXXXXXXXXXXXXXXXXXX1)]==>Z   (27) correct
       1 0 0 A+S[USER;cXXXXXXXXXXXXXXXXXXXXXXXXXX1:(...)]==>Z                                (28) FATAL #999
   (if 27-character species names were intended, the limit is one short). The converter now keeps
   species names to 26 characters, so 5,228 species in 422 packs got their SBML names back.
   Packs built this way need the 2026-10-26 15:30 build or later.

STILL PRESENT
2. LONG USER EXPRESSIONS ARE MIS-READ - FATAL #999 or a silently WRONG value. This is the cause of every
   FATAL #999 seen in the collection (35 packs) and of many stalls. Deck: A=1 free, B=0.5 pinned; 1 s;
       1  0  0  A+S[USER;cA,cB:(EXPR)]==>Z
   (a) EXPR = 16 copies of 0.001*cB^(cB+0.5) joined by '+' (287 characters). Correct rate 0.008 -> A(1)=0.992032.
       Kintecus gives A(1)=0.993024 = exp(-0.007): two terms are dropped, no message. 18, 20, 24, 30, 36 copies
       likewise lose exactly two terms. 14 copies (251 characters) -> FATAL #999. 10 copies (179) correct.
   (b) Any '^' in the expression: correct up to 248 characters, FATAL #999 or wrong from ~252 on
       (0.001*cB^2 followed by +0.0001 terms padded to length L: L=248 right, 252 FATAL, 256-268 wrong).
       With more '^' the limit drops a little (8 powers: 244 already fails).
   (c) Without '^' long expressions work (constant/variable terms right at 740 characters), but some
       lengths still give FATAL #999 depending only on where the characters fall: the same 573-character
       expression with 0..9 extra trailing zeros on one constant (same value) -> OK OK 999 999 OK OK 999 OK OK OK.
       The case reported earlier (BIOMD0000000013: a 400-character expression evaluating to 0; changing
       0.00064 -> 0.0006 fixed it) is the same issue.
   It looks like a fixed ~255-character work buffer in the expression handling ('^' path), plus an
   off-by-one/boundary effect elsewhere. Converter workaround: see README (shorter text, '^'-free powers,
   short names, per-line evaluation check and re-padding); a longer buffer in Kintecus would remove it.
3. At most 3 reactant entries per reaction line: A+B+C+D==>E (plain or with S[USER...]) gives FATAL #39
   "exceeded the number of reactants". (Not line length - an earlier note blaming long lines was wrong.)
   9 packs hit this; the converter now keeps lines to 3 reactants.
6. A species referenced only inside a USER expression (in no reaction line) stalls the run or gives
   FATAL #999 unless it is pinned with SPECIES column 7.
7. A USER rate with a time discontinuity (step built from abs()) stalls the integrator at the step;
   Constant File profiles handle the same discontinuity cleanly.
8. Negative zero-order USER rates crawl/stall; negative product coefficients (B+-1A) are not parsed.
9. '0 NULL+S[USER...]==>X' whose rate depends on X hides that dependence from the Jacobian; on a stiff
   system (Oregonator) the result is wrong. Writing X on both sides (gross stoichiometry) fixes it.

OPEN - not diagnosed
10a. 2026-10-30: the Cash-Karp integrator (-INT:2, no Jacobian) was added as a third attempt on the 170 packs that
    ran out of time with -INT:1 and -INT:5 (300 s; 600 s without the projected-time cut-off for the 19 closest).
    It finishes some of them quickly where the Jacobian-based integrators crawl (e.g. BIOMD0000000496/497 in
    2 s); on strongly stiff packs it stalls, or finishes with NaN/Inf (e.g. BIOMD0000000214, 571). All 170 packs
    that ran out of time use USER rate lines (none of the 225 pure mass-action packs did), consistent with the
    interpreted USER evaluation, multiplied by finite-difference Jacobians, being the cost.
10b. FATAL #95 with -INT:2 (no message) on BIOMD0000000127, 582, 801; the packs are the repro decks.
10c. 'SOFTWARE ERROR AT UFHS88' right after 'Running Simulation' on 13 large USER packs (470-2,691 lines,
    400+ USER lines; e.g. BIOMD0000000332, 603-610) - likely an internal table limit.
10. Integrator behavior on the folder-2 packs (deck verified, Kintecus too slow): all 200 were rerun with
    nothing else running. -INT:1 does advance but would need 2,000-26,000 s (300 s allowed); both
    -INT:5:59:1e-20 and plain -INT:5 produce no output row at all within 30 s on almost every one
    (BIOMD0000000018: none in 140 s with -INT:5 at accuracy 1e-8 or 1e-10, while -INT:1 reached t=500).
    Packs whose deck is verified correct but where Kintecus stalls at t=0 or is far too slow with both
    -INT:1 and -INT:5:59:1e-20 (folder 2; FINAL_RESULTS_TABLE notes give the attempt details), and packs
    where Kintecus completes but differs although the deck is verified (also folder 2). Example to start
    with: BIOMD0000000202 (Ca2+ / ER channel; stalls at t=0 on both integrators).

INTERFACE NOTES
11. With -i, a FATAL number is written only to KERR.OUT (KDONE.OUT: "FATAL ERROR ENCOUNTERED!!").
    KRUN.TMP only says "Kintecus is running...." - a current simulated time there would let external tools
    tell "slow" from "stuck" exactly (CONC.TXT growth is used now). -i needs no key presses at all.
12. Kintecus enforces at least ~61 output rows; Constant File / -SET / -ADD changes take effect at the
    first output step at or after the profile time, so event times must lie on the sampling grid to be
    exact - and the grid is reached by accumulating the sampling interval, so an interval like 1/7 written
    exactly (0.14285714285714285) makes 7 steps fall a hair short of t = 1 and the change lands one step
    late (BIOMD0000000818: 0.002% -> 99.9% off). The converter now writes the interval rounded UP in the
    10th digit (0.1428571429). A tolerance of ~1e-12 x interval in Kintecus' "has the profile time been
    reached" test would remove the trap. Accumulation does it too: 100 x 0.1 = 9.99999999999998 < 10, so
    with any interval that is not exact in binary a change at t=10 also lands one step late. The converter
    now writes the interval a hair long when a deck has profiles (0.1 -> 0.1000000002); this moved
    BIOMD0000000001, 120, 318, 404, 527, 563, 601, 806, 818, 820, 918, 955 and 960 to PASS.
15. "Run once and hold" profiles (negative spacing) with -SET or -ADD keep re-applying the LAST entry every
    period after the profile ends: a -SET profile 0, 5 at spacing -2 sets 5 at t = 2, 4, 6, 8, 10, ...;
    -ADD keeps adding. Plain (pinned) profiles hold correctly. The converter ends -ADD profiles with a 0;
    -SET cannot be worked around in the deck. BIOMD0000000077, 281, 337, 1019 and 1020 differ for this
    reason (their decks reproduce the SBML; Kintecus diverges at twice the profile spacing).
    Since 2026-11-05 the packs use Time/value-pair files (item 21), which apply each value once; the five
    -SET packs affected (0077, 0281, 0337, 1019, 1020) now match the SBML within 0.002%.
16. (2026-11-01) A variable driven below zero stalls the integration. Repro: repros/NEG_FROM_ZERO_stall/
    (A = 0, B = 1.25; -0.8275 source==>A; 0.678 B==>sink; 1 source+S[USER;cB:(0.662*cB)]==>A; exact
    A(10) = -7.0559). -INT:1, -INT:2 and -INT:5 produce no row after t = 0. -INT:5:10:0 (DASPK without the
    non-negativity constraint) gives A(10) = -7.05589, then loops at the end ("AFTER DASPK, IDID= -33",
    TEMP 9.99999999999999 vs RTEMP 10) and never writes KDONE.OUT. Many SBML models have state variables
    that go to or below zero (parameters with rate rules such as deviations from a basal level); the packs
    now write those terms exactly (signed lines) and also try -INT:5:10:0. See repros/README_NEGATIVE_VALUES.txt.
17. (2026-11-01) -INT:5 and -INT:2 can run to the end time with NaN in CONC.TXT and only Warning #6 in
    KWARN.OUT (e.g. BIOMD0000000180, 215, 568, 675, 896, 961 before the 11-01 converter changes). A FATAL
    error there would be clearer.
18. (2026-11-01) OUTC[] with a one-character factor at the END of a term fails while reading SPECIES.DAT:
    OUTC[As;-source*2] or OUTC[As;-source*5] -> "*ERROR* Soft error u82n (M+1).GE.FIRSTN(temp2)" and no
    output; OUTC[As;-source*20], OUTC[As;-2*source] and OUTC[As;-source*2.0] work. The packs write the factor
    first (OUTC[X_o;-2*source]).
19. (2026-11-02) A Constant File profile can hold at most 20,002 values (20,003 lines with the spacing line):
    one more value gives FATAL Error #94 while reading SPECIES.DAT. The packs keep profiles at or under
    20,001 values and stop a profile just after its last event (negative spacing holds the last value).
    BIOMD0000001029 needs 36,002 values, so it stays blocked until the limit is raised (40,000 would cover it).
    Since 2026-11-05 BIOMD0000001029 uses a Time/value-pair file instead and matches the SBML (0.00047%).
20. (2026-11-02) Far into a long profile, an -ADD entry can be applied one interval late (profile clock in
    limited precision?). Repro: repros/PROFILE_ADD_TIMING (spacing -0.1, 1.0 every 500th entry): the steps
    show at 100.2-300.2 and 1250.2-1650.2 instead of x.1. BIOMD0000000537: the 1344 h dose lands 0.1 h late
    (99% difference on one row; 0.6% elsewhere), so it is reported as MISMATCH.
    Since 2026-11-05 the packs use Time/value-pair files (item 21); BIOMD0000000537 now matches (0.00046%).
21. FIXED in the build compiled 2026-11-05 08:40:46 (md5 9b45994d...): Time/value-pair Constant Files (first
    line contains "Time", pairs until END). On the 11-03 builds only the FIRST pair was applied, a single pair
    could be skipped (14 s pair lost in a 60 s run) and a pair below the first step (1e-12) was never applied.
    Now correct under -INT:1, -INT:5 and -INT:2 for -ADD, -SET, 1,000 pairs in one file and 50 pair files.
    Still open (minor): a PINNED pair file (no -ADD/-SET) reads -987654000 before its first pair instead of
    the SPECIES.DAT value; the packs never use pinned pair files. Repros: repros/TIMEPAIRS_*,
    repros/README_1105_BUILD.txt.
23. (2026-11-05, question) After each Time/value pair Kintecus writes a row at the pair time and restarts the
    output grid from it (pair at 14, step 0.98: rows 13.77, 14.0, 14.98, ... instead of 13.77, 14.75, ...);
    spacing-format profiles keep the grid. Values are correct; row times move. The references in this
    release are written at the times Kintecus prints (if the grid changes, kintcheck interpolates between
    rows, which is inexact across a dose). Repro: repros/TIMEPAIRS_OUTPUT_GRID (pairs/ and spacing/).
24. (2026-11-07) -INT:1 never applies Time/value-pair changes in one deck: BIOMD0000000554 (3 switch files at
    t=200, 500, 723). Under -INT:1 the switch species keep their initial values all run and the output jumps
    197.5 -> 212.077 (no row at the pair time); -INT:5 applies them and matches the SBML (0.027%). Small decks
    with the same pair files switch correctly. Repro: repros/PAIR_SKIP_INT1_BIOMD0000000554 (README.txt).
25. (2026-11-07) Kintecus's stiff integrators are much slower than CVODE/LSODA on a few correct decks:
    BIOMD0000000506 (56 species, 108k characters of USER: LSODA 3 s, CVODE 0.02 s, -INT:1 reaches t=220 of 1000
    in 150 s regardless of sampling or maximum step), BIOMD0000000608 (LSODA 0.4 s; -INT:1 t=41,000 of 100,000
    in 600 s, view.txt "Stepsize underflow in stiff sub" x1139), BIOMD0000000603 after t=5425 (-INT:1 stops where
    the cell volume sits on the turgor threshold, a kink; -INT:5 passes it). Long USER lines cost time per
    evaluation: 160 lines of 366 characters added 10 s to a 0.8 s run (slowdiag/jac/LONG*).
22. FIXED in the build compiled 2026-11-03 16:34:03: USER rates were re-parsed on every call on the 11-03
    12:39/13:35 builds (interpreter cache key 1,024 characters < 2,159-character caller buffer; earlier builds
    thrashed a 64-entry round-robin cache). New hashed cache: same output bit-for-bit, 33-64x faster on USER
    models (BIOMD0000000010, 0100, 0255, 1028). Repro: repros/SPEED_USER_BIOMD0000000010.
13. Under Wine, Kintecus.exe is re-parented away from its launcher; killing the launcher leaves it running
    (the harness here now kills it by working directory). Not an issue on Windows.
14. FIXED in the build compiled 2026-11-03 16:34:03: the reaction field was limited to 897 characters; it is
    now 2,153 (2,154 -> FATAL #24). Repros: repros/REACTION_FIELD_897, REACTION_FIELD_898. The packs keep
    reaction comments on a '#' line above each reaction (longest reaction field in the collection: 2,148).
"""

# DESCRIPTIONS.txt
out = io.StringIO()
out.write(f"Pack descriptions - {len(rows)} curated BioModels entries ({today})\n")
out.write("Folder meanings and scoring: see README.txt. The same notes are in FINAL_RESULTS_TABLE.csv.\n\n")
for r in rows:
    R = res[r["id"]]
    out.write("=" * 100 + "\n")
    out.write(f"{r['id']}  {r['name']}\n")
    out.write(f"  pack: {r['kintpack'] or '(none)'}   folder: {r['folder']}   SUCCESS: {r['success_under_4pct']}   within 0.1%: {r['within_0.1pct']}\n")
    out.write(f"  Kintecus vs SBML: {r['kintecus_match_error_percent'] or 'n/a'}%   deck vs SBML: "
              f"{r['deck_match_error_percent'] or 'n/a'}%   (max |diff| / species peak)   status {r['status']}\n")
    if r["kintpack"]:
        out.write(f"  run: {R.get('switches', '')}   0..{r['duration']} ({r['duration_source']}), sampling "
                  f"{R.get('sampling', 0):.6g}; {r['sbml_reactions']} SBML reactions -> {r['kintecus_lines']} "
                  f"Kintecus lines, {r['species']} species\n")
    if r["kintecus_runs"]:
        out.write(f"  Kintecus runs: {r['kintecus_runs']}\n")
    out.write("  note: " + r["notes"] + "\n")
    nc = [t for s, t in R.get("flags", []) if s == "blocker"]
    if nc:
        out.write("  NOT converted:\n")
        for t in nc[:12]:
            out.write(f"    - {t}\n")
        if len(nc) > 12:
            out.write(f"    ... and {len(nc) - 12} more (results_detail.csv)\n")
open(f"{B}/DESCRIPTIONS.txt", "w").write(out.getvalue())

# results_detail.csv
detail = []
for r in rows:
    R = res[r["id"]]
    blockers = [t for s, t in R.get("flags", []) if s == "blocker"]
    detail.append({**r, "run_switches": R.get("switches", ""), "sampling": R.get("sampling"),
                   "user_lines": R.get("n_user_lines"), "long_user_lines": R.get("long_user_lines", ""),
                   "long_user_lines_misread_then_fixed": R.get("long_user_lines_repadded", ""),
                   "not_converted": " | ".join(blockers)[:2000],
                   "conversion_notes": " | ".join(f"[{s}] {t}" for s, t in R.get("flags", []) if s != "blocker")[:3000]})
with open(f"{B}/results_detail.csv", "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(detail[0])); w.writeheader(); w.writerows(detail)

bat = open(f"{B}/kintcheck/RUN_ALL_CHECKS.bat").read().replace("\r\n", "\n").replace("\n", "\r\n")
_mb = lambda paths: sum(os.path.getsize(p) for p in paths) / 2 ** 20
_fold = {f: [f"{B}/packs/{r['kintpack']}" for r in packs if r["folder"] == f] for f in
         ("1_success_under_4pct", "2_deck_verified_kintecus_too_slow_or_differs", "3_off_partial_or_unverified")}
_ref_mb = _mb([f"{B}/references/{r['id']}.csv" for r in packs if r["id"] in refs])
CONTENTS = f"""CONTENTS - biomodels_kintpacks_final.zip   ({today})
================================================================================================
Start with README.txt, then FINAL_RESULTS_TABLE.csv.

PACKS ({len(packs)} .kintpack files, sorted by result)
  1_success_under_4pct/                          {len(_fold['1_success_under_4pct']):4d} packs  {_mb(_fold['1_success_under_4pct']):6.1f} MB
      Kintecus output matches the original SBML to under 4% (most within 0.1%; see within_0.1pct).
  2_deck_verified_kintecus_too_slow_or_differs/  {len(_fold['2_deck_verified_kintecus_too_slow_or_differs']):4d} packs  {_mb(_fold['2_deck_verified_kintecus_too_slow_or_differs']):6.1f} MB
      The deck reproduces the SBML when run on its own, but Kintecus is too slow, stops, or gives different numbers.
  3_off_partial_or_unverified/                   {len(_fold['3_off_partial_or_unverified']):4d} packs  {_mb(_fold['3_off_partial_or_unverified']):6.1f} MB
      Off by 4% or more, has parts that couldn't be converted, or has no SBML result to check against.
  The other {cnt['(no pack)']} of the {len(rows)} models produced no pack; FINAL_RESULTS_TABLE says why.

RESULTS AND DOCUMENTATION
  README.txt                What's in the zip, how to load a pack, how the match is scored, how to run the
                            self-check, how to convert your own SBML files, and the conversion conventions.
  CONTENTS.txt              This file.
  FINAL_RESULTS_TABLE.csv   One row per model ({len(rows)} rows), opens in Excel. Columns: success yes/no at 4%, within 0.1% yes/no,
                            Kintecus match %, deck match %, status, folder, a note on why it failed or how far
                            off it is, what each integrator attempt did, number of unconverted parts, duration,
                            and whether the pack was rebuilt with the latest fixes.
  FINAL_RESULTS_TABLE.txt   The same table as readable text, with a summary of counts at the top.
  DESCRIPTIONS.txt          One block per model: results, switches, run length, sizes, notes, and the full list
                            of anything not converted.
  KINTECUS_ISSUES.txt       Kintecus V2027 problems found, each with a minimal repro: fixed in your build, still
                            present (long USER expressions behind #999, the 3-reactant limit behind #39, function
                            names inside variable names, ...), open, and interface notes.
  results_detail.csv        The full data behind the table, including conversion notes and the not-converted
                            list per model.

SELF-CHECK ON YOUR PC
  RUN_ALL_CHECKS.bat        Double-click to run every pack through your Kintecus.exe and compare each against its
                            SBML reference. Set the Kintecus.exe path on one line, or answer the prompt.
  RUN_SLOW_JOBS.bat         Runs the packs in SLOW_JOBS.csv (those that ran out of time here) with no time limit,
                            each with its recommended integrator, and checks them (slow_jobs_report.csv).
  SLOW_JOBS.csv             The slow packs: order, recommended switches, run-time estimate, deck verified,
                            other integrators to try, what each attempt did here.
  SLOW_JOBS_README.txt      What to do on your PC to run and finish those jobs.
  NOT_TIME_PROBLEMS.csv     Packs that stopped on a Kintecus error (FATAL #999/#167/#95, UFHS88), not on time.
  kintcheck.py              The checker the .bat runs (Python 3.8+, no extra packages). Writes
                            kintcheck_report.csv / .txt, including the result each pack gave when it was built.
  references/               {n_ref} CSV files of reference output from the original SBML ({_ref_mb:.0f} MB, the bulk of
                            the zip). The other {len(packs) - n_ref} packs have no reference (their SBML can't be simulated by
                            libRoadRunner, or nothing in it changes); the checker reports those as NO_REFERENCE.

CONVERT YOUR OWN SBML FILES
  sbml2kintpack/sbml2kintpack.py   One SBML file in -> <name>.kintpack + <name>_deck\\ (MODEL.DAT, SPECIES.DAT,
                                   PARM.DAT, profile files). Duration from --t-end, a .sedml next to the model,
                                   or the dynamics. --kintecus PATH checks long USER lines in Kintecus itself.
                                   Needs: pip install python-libsbml sympy (libroadrunner numpy recommended).
  sbml2kintpack/sbml2kin_batch.py  The converter it uses (the same one that built this collection).
  sbml2kintpack/release_settings.csv  Duration, sampling, switches, accuracy, deck variant of every released pack.
  sbml2kintpack/MAKE_KINTPACK.bat  Drag SBML files onto it (set the Kintecus.exe path at the top for the check).

RUN A PACK FROM THE COMMAND LINE
  runkintpack/runkintpack.py       Unpacks any .kintpack (all sheets and files, tab-delimited) into a run folder
                                   and runs Kintecus.exe there; the .kintpack is never changed. -k switches,
                                   -p Kintecus.exe, -i international inputs, --set sheet.field=value overrides
                                   (species, parm, model lines, therm, transport, fitdata/profile files) with
                                   checking, --list, --check, -u, --batch. Python 3.8+, no extra packages.
  runkintpack/runkintpack.bat      runkintpack.bat PACK.kintpack [options]  (finds Python for you)
  runkintpack/RUNKINTPACK_MANUAL.txt  the manual: options, the --set syntax, examples, messages

EXPORT TO SBML
  kintpack2sbml/kintpack2sbml.py   .kintpack or Kintecus deck -> SBML L3V2 (+ SED-ML, + report of everything SBML
                                   cannot hold; --format xml/both for a lossless Kintecus XML). Needs python-libsbml.
  kintpack2sbml/kintpack2sbml.bat  kintpack2sbml.bat INPUT [options]
  kintpack2sbml/KINTPACK2SBML_MANUAL.txt  the manual: options, how Kintecus maps onto SBML, examples, test results

converter/  (reference: the Linux scripts that built and verified the collection; they need the BioModels
            mirror, libRoadRunner and Wine, so they are not meant to be run on your PC)
  sbml2kin_batch.py       the SBML -> Kintecus converter
  run_model.py, driver.py convert + run + verify one model / a batch
  deck_sim.py             runs a deck on its own (SciPy, Kintecus semantics) - the independent deck check
  validate_packs.py       checks each pack against the deck that was run
  gen_refs.py             builds the reference files
  make_final_table.py     builds FINAL_RESULTS_TABLE
  package_v3.py           builds this zip
  probe.py                the Kintecus run harness (Wine)
"""
open(f"{B}/CONTENTS.txt", "w").write(CONTENTS)

zpath = f"{B}/biomodels_kintpacks_final.zip"
with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
    z.writestr("README.txt", README.replace("\n", "\r\n"))
    z.writestr("KINTECUS_ISSUES.txt", ISSUES.replace("\n", "\r\n"))
    z.writestr("DESCRIPTIONS.txt", out.getvalue().replace("\n", "\r\n"))
    z.write(f"{B}/FINAL_RESULTS_TABLE.csv", "FINAL_RESULTS_TABLE.csv")
    z.writestr("FINAL_RESULTS_TABLE.txt", open(f"{B}/FINAL_RESULTS_TABLE.txt", encoding="utf-8").read().replace("\n", "\r\n"))
    z.write(f"{B}/results_detail.csv", "results_detail.csv")
    z.writestr("RUN_ALL_CHECKS.bat", bat)
    z.writestr("CONTENTS.txt", CONTENTS.replace("\n", "\r\n"))
    z.write(f"{B}/sbml2kintpack/sbml2kintpack.py", "sbml2kintpack/sbml2kintpack.py")
    z.write(f"{B}/sbml2kintpack/release_settings.csv", "sbml2kintpack/release_settings.csv")
    z.write(f"{B}/runkintpack/runkintpack.py", "runkintpack/runkintpack.py")
    z.write(f"{B}/kintpack2sbml/kintpack2sbml.py", "kintpack2sbml/kintpack2sbml.py")
    for fn in ("kintpack2sbml.bat", "KINTPACK2SBML_MANUAL.txt"):
        z.writestr(f"kintpack2sbml/{fn}",
                   open(f"{B}/kintpack2sbml/{fn}").read().replace("\r\n", "\n").replace("\n", "\r\n"))
    z.writestr("runkintpack/RUNKINTPACK_MANUAL.txt",
               open(f"{B}/runkintpack/RUNKINTPACK_MANUAL.txt").read().replace("\r\n", "\n").replace("\n", "\r\n"))
    z.writestr("runkintpack/runkintpack.bat",
               open(f"{B}/runkintpack/runkintpack.bat").read().replace("\r\n", "\n").replace("\n", "\r\n"))
    z.write(f"{B}/sbml2kin_batch.py", "sbml2kintpack/sbml2kin_batch.py")
    z.writestr("sbml2kintpack/MAKE_KINTPACK.bat",
               open(f"{B}/sbml2kintpack/MAKE_KINTPACK.bat").read().replace("\r\n", "\n").replace("\n", "\r\n"))
    z.write(f"{B}/kintcheck/kintcheck.py", "kintcheck.py")
    for fn in ("MODEL.DAT", "SPECIES.DAT", "PARM.DAT"):
        z.write(f"{B}/repro_neg/NEG_FROM_ZERO_stall/{fn}", f"repros/NEG_FROM_ZERO_stall/{fn}")
    for d in ("PROFILE_20002_VALUES", "PROFILE_20003_VALUES"):
        for fn in ("MODEL.DAT", "SPECIES.DAT", "PARM.DAT", "PROF.TXT"):
            z.write(f"{B}/repro_null/{d}/{fn}", f"repros/{d}/{fn}")
    for fn in ("MODEL.DAT", "SPECIES.DAT", "PARM.DAT", "ADD.TXT"):
        z.write(f"{B}/repro_null/PROFILE_ADD_TIMING/{fn}", f"repros/PROFILE_ADD_TIMING/{fn}")
    for d in ("REACTION_FIELD_897", "REACTION_FIELD_898", "SPEED_USER_BIOMD0000000010", "TIMEPAIRS_ADD_3PAIRS",
              "TIMEPAIRS_PINNED_2PAIRS", "TIMEPAIRS_RUN_60S", "TIMEPAIRS_RUN_70S", "TIMEPAIRS_BIOMD0000000439",
              "TIMEPAIRS_OUTPUT_GRID/pairs", "TIMEPAIRS_OUTPUT_GRID/spacing", "PAIR_SKIP_INT1_BIOMD0000000554"):
        for fn in sorted(os.listdir(f"{B}/repro_null/{d}")):
            if os.path.isfile(f"{B}/repro_null/{d}/{fn}") and fn.upper() not in ("KINTECUS.EXE",):
                z.write(f"{B}/repro_null/{d}/{fn}", f"repros/{d}/{fn}")
    z.writestr("repros/README_1105_BUILD.txt", open(f"{B}/repro_null/README_1105_BUILD.txt").read().replace("\r\n", "\n").replace("\n", "\r\n"))
    z.writestr("repros/README_1103_BUILD.txt", open(f"{B}/repro_null/README_1103_BUILD.txt").read().replace("\r\n", "\n").replace("\n", "\r\n"))
    z.writestr("repros/README_REPROS.txt", open(f"{B}/repro_null/README_REPROS.txt").read().replace("\r\n", "\n").replace("\n", "\r\n"))
    z.writestr("repros/README_NEGATIVE_VALUES.txt",
               open(f"{B}/repro_neg/README_NEGATIVE_VALUES.txt").read().replace("\n", "\r\n"))
    z.writestr("RUN_SLOW_JOBS.bat", open(f"{B}/kintcheck/RUN_SLOW_JOBS.bat", newline="").read())
    z.writestr("RUN_REGRESSION.bat", open(f"{B}/kintcheck/RUN_REGRESSION.bat", newline="").read())
    z.write(f"{B}/kintcheck/SLOW_JOBS.csv", "SLOW_JOBS.csv")
    z.write(f"{B}/kintcheck/NOT_TIME_PROBLEMS.csv", "NOT_TIME_PROBLEMS.csv")
    z.writestr("SLOW_JOBS_README.txt", open(f"{B}/kintcheck/SLOW_JOBS_README.txt").read().replace("\n", "\r\n"))
    for r in packs:
        if r["id"] in refs:
            z.write(f"{B}/references/{r['id']}.csv", f"references/{r['id']}.csv")
    for fn in ("sbml2kin_batch.py", "run_model.py", "deck_sim.py", "driver.py", "validate_packs.py", "gen_refs.py",
               "make_final_table.py", "package_v3.py"):
        z.write(f"{B}/{fn}", f"converter/{fn}")
    z.write("/home/claude/kin/probe/probe.py", "converter/probe.py")
    for r in packs:
        z.write(f"{B}/packs/{r['kintpack']}", f"{r['folder']}/{r['kintpack']}")
print(README[:900])
print("zip MB:", round(os.path.getsize(zpath) / 2 ** 20, 1), " packs:", len(packs), " refs:", n_ref, " success:", n_ok)
