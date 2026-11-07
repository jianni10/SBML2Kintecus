BioModels -> Kintecus .kintpack collection
==========================================
Built 2026-11-07. Source: the curated BioModels set as mirrored in github.com/sys-bio/temp-biomodels
(entries up to BIOMD0000001080, mirror dated June 2025).

1072 curated models -> 1045 .kintpack files.
SUCCESS (Kintecus output under 4% of the original SBML; rule set 2026-10-30): 943 packs
  (922 of them within 0.1%, the earlier rule; column within_0.1pct in the table).
Packs shown to reproduce the SBML (SUCCESS, or folder 2 where an independent SciPy run of the deck matches
  within 0.1% but Kintecus is too slow or differs): 947.
1072 entries (including every former #39/#999 pack) were rebuilt with the latest converter
fixes; the others keep their earlier build and result (column built_with_latest_fixes in the table).

START HERE: FINAL_RESULTS_TABLE.csv (open in Excel) or FINAL_RESULTS_TABLE.txt - one row per model with the
match %, SUCCESS yes/no, and a note on why it failed or how far off it is.

WHAT IS IN THIS ZIP   (full table of contents: CONTENTS.txt)
  1_success_under_4pct/                           943  Kintecus output matches the SBML to under 4%
  2_deck_verified_kintecus_too_slow_or_differs/     4  the deck reproduces the SBML (independent
                                                       check), but Kintecus is too slow, stops, or differs
  3_off_partial_or_unverified/                     98  off by 4% or more, too slow, parts not convertible,
                                                       or no SBML reference to check against
  (no pack)                                        27  nothing convertible / converter error
  FINAL_RESULTS_TABLE.csv / .txt   the results table (see above)
  RUN_ALL_CHECKS.bat               re-run every pack on YOUR PC and compare with the SBML reference (below)
  RUN_SLOW_JOBS.bat                run the packs that did not finish here, with no time limit (SLOW_JOBS_README.txt)
  RUN_REGRESSION.bat               after a Kintecus change: re-run the passing packs and list the ones that now fail
                                   (regression_report_FAILED.csv; FAILED = over 4%, WORSE = error grew > 10x)
  kintcheck.py                     the checker the .bat runs (plain Python 3, no packages needed)
  (original SBML + SED-ML inputs: separate download biomodels_original_sbml.zip, one folder per entry)
  references/                      SBML reference output for 1029 packs (libRoadRunner), used by kintcheck.py
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
        if "%KINTECUS_EXE%"=="" set "KINTECUS_EXE=C:\Kintecus\Kintecus.exe"
     (or put Kintecus.exe next to the .bat, or answer the prompt when it asks).
  4. Double-click RUN_ALL_CHECKS.bat. For each pack it writes the DAT/profile files into
     kintcheck_work\<ID>\ exactly as the GUI does, runs Kintecus with the pack's switches plus -i, watches
     KDONE.OUT / KERR.OUT / CONC.TXT, and compares CONC.TXT with references\<ID>.csv.
  5. Results: kintcheck_report.txt (summary) and kintcheck_report.csv (per pack: PASS/OFF/STALLED/...,
     match error %, worst species, Kintecus message, and the result this pack gave when it was built, so a
     new Kintecus version's changes stand out). Failed runs stay in kintcheck_work\ for inspection.
  Options (add after the .bat name in a command prompt, or edit the .bat):
     --folder 1_success          only packs in matching folders     --only BIOMD0000000010 ...  only these IDs
     --stall 120 --cap 1800      more patience (seconds)            --tol 0.04                   threshold (4%)
     --switches "-ig:mass -INT:5:59:1e-20 -obeymaxint"  force other switches     --show  show Kintecus windows
  A full run of all packs takes hours (each stuck pack waits --stall seconds; slow ones up to --cap).
  Start with  RUN_ALL_CHECKS.bat --folder 1_success  (about 943 packs).

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
  Drag SBML files onto sbml2kintpack\MAKE_KINTPACK.bat, or from a command prompt:
     python sbml2kintpack.py model.xml                         duration from a .sedml next to the model
     python sbml2kintpack.py model.xml --t-end 500 --points 1000
     python sbml2kintpack.py model.xml --kintecus "C:\Kintecus\Kintecus.exe"
  Output next to the model: <name>.kintpack (load in the GUI) and <name>_deck\ (MODEL.DAT, SPECIES.DAT,
  PARM.DAT, profile files). --kintecus runs the long-USER-expression check in Kintecus itself and repairs
  lines Kintecus mis-reads (recommended; MAKE_KINTPACK.bat does it when it finds Kintecus.exe - set the
  path at the top of the .bat). For a BioModels entry (ID in the file or folder name, or the model's own
  annotation) the duration, sampling, switches, accuracy and deck variant (shifted states etc.) of the released
  pack are taken from sbml2kintpackelease_settings.csv: converting any of the 890 SUCCESS models from
  biomodels_original_sbml.zip this way was checked here to give the verified deck (838 identical; 52 written by
  the current converter instead of an older one, all 52 run in Kintecus and PASS). --no-release-settings skips
  the table. States that go below zero are written shifted with an OUTC[] output (--no-shift: signed lines).
  Other SBML files are converted the same way but are not verified - check the description's conversion notes
  and "NOT converted" list.

RUNNING A PACK FROM THE COMMAND LINE (runkintpack/)
  runkintpack.bat PACK.kintpack [options]        (or: python runkintpack.py ...; Python 3.8+, no packages)
  Writes every sheet and file of the pack - MODEL.DAT, SPECIES.DAT, PARM.DAT, THERM.DAT, TRANSPORT.DAT,
  fitdata and perturbation/profile files - tab-delimited into <pack>_run\ (or -o DIR) and runs Kintecus.exe
  there. The .kintpack itself is never changed.
     -k "SWITCHES"   override the pack's Kintecus switches       -p PATH   Kintecus.exe (overrides the pack's)
     -i              international inputs: 4,34 -> 4.34 (--set fields are then separated by ';')
     --set species.H2O.initialconcentration=4.34      --set parm.temperature=897
     --set model.12="3.4e10,0,0,A+B==>C,my comment"  --set model.12.A=3.4e10
     --set fitdata001.5="0,1.25,3.5"                  --set therm.3.2=NASA     --set transport.4="..."
     --list (show settable fields)  --check (validate --set only)  -u (unpack only)
     --batch [--stall S] [--cap S]  (unattended: hidden console, stops stuck runs)
  runkintpack.bat --help lists everything; runkintpack\RUNKINTPACK_MANUAL.txt has the full manual with examples.

EXPORTING TO SBML (kintpack2sbml/)
  kintpack2sbml.bat INPUT [options]   INPUT = a .kintpack, a Kintecus deck folder, or its MODEL.DAT
  Writes NAME.xml (SBML L3V2), NAME.sedml (time course) and NAME_kintecus_unsupported.txt (everything SBML
  cannot hold, with the original text; also in an XML comment and the model notes). --format both adds a
  lossless NAME.kintecus.xml. --validate runs the libSBML check. Needs: pip install python-libsbml.
  Full manual: kintpack2sbml\KINTPACK2SBML_MANUAL.txt

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
