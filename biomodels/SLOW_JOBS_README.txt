SLOW JOBS - BioModels packs that did not finish within the time limits used here, and how to finish them
=====================================================================================================  (2026-11-13)

WHICH PACKS
  SLOW_JOBS.csv lists 11 packs. Each one ran out of time (too slow, or no new output row for the stall
  window) on -INT:1 and -INT:5, 300 s each. -INT:2 (explicit Cash-Karp Runge-Kutta) is no longer used: it
  cannot integrate stiff systems (runs end with NaN/Inf), so only -INT:1 and -INT:5 are suggested. Kept run
  folders have VIEW.TXT / KWARN.OUT over 20 MB cut to their first and last 1,000 lines (*_head_tail.txt).
    - 1 have a verified deck: an independent SciPy run of the same pack matches the SBML, so the pack is right
      and only Kintecus's run time stands in the way. These come first in the list.
    - the rest are not confirmed (the independent run was also off or could not finish), so a finished run may
      still differ from the SBML.
  Packs that failed on a Kintecus error rather than on time (FATAL #999, FATAL #167, FATAL #95, SOFTWARE ERROR AT
  UFHS88) are in NOT_TIME_PROBLEMS.csv; more time will not help them.

WHAT TO DO ON YOUR PC
  1. Unzip biomodels_kintpacks_final.zip into a folder, e.g. C:\BioModels (keep the folder structure: the packs in
     their 1_/2_/3_ folders, references\, kintcheck.py, RUN_SLOW_JOBS.bat and SLOW_JOBS.csv side by side).
  2. Python 3.8 or newer must be installed (python.org); no extra packages.
  3. Open RUN_SLOW_JOBS.bat in Notepad and set the Kintecus.exe path on the line
         if "%KINTECUS_EXE%"=="" set "KINTECUS_EXE=C:\Kintecus\Kintecus.exe"
     (use the 2026-11-05 08:40:46 build or later: earlier builds misapply Time/value-pair Constant Files and
     fail on packs with reaction lines over 897 characters or more than 52 variables on a USER line).
  4. Double-click RUN_SLOW_JOBS.bat. It runs the jobs one after another, fastest estimate first, each with the
     switches in SLOW_JOBS.csv, with NO time limit and NO stall limit, prints progress every minute, and compares
     each finished run with the SBML reference (PASS = under 4% of peak; --tol 0.001 for the strict 0.1% check).
  5. Results go to slow_jobs_report.txt / .csv; every run folder is kept in kintcheck_work\<ID>\.
  6. Send me slow_jobs_report.csv and I will fold the results into FINAL_RESULTS_TABLE, the plots and the Q&A.

  The whole list can take days. To run part of it:
     RUN_SLOW_JOBS.bat --only BIOMD0000000316 BIOMD0000000401      just these packs
     RUN_SLOW_JOBS.bat --cap 14400                                 stop any single job after 4 h
                       (with --cap, a job is also stopped once its projected time passes the cap;
                        add --no-projection to give every job the full cap)
     RUN_SLOW_JOBS.bat --stall 1800                                stop a job after 30 min without a new row
  Several copies can run at the same time on a multi-core PC if each gets its own --only list and --name, e.g.
     RUN_SLOW_JOBS.bat --only BIOMD0000000017 --name part1
     RUN_SLOW_JOBS.bat --only BIOMD0000000049 --name part2

  One integrator for every job (e.g. a final DASPK check), each job given the full 30 minutes:
     RUN_SLOW_JOBS.bat --switches "-ig:mass -INT:5:10:0 -obeymaxint" --cap 1800 --no-projection --name slow_jobs_int5
  (-INT:5:10:0 = DASPK without its non-negativity constraint: packs written with signed lines let a variable go
   below zero exactly where the SBML does, and plain -INT:5 stalls there; use "-ig:mass -INT:5 -obeymaxint" for
   plain DASPK.)

  To run one by hand instead: load the pack in the Kintecus GUI ("Load Package..."), put the switches from the
  'switches' column in the switches box, and run; or from a command prompt:
     runkintpack\runkintpack.bat 2_deck_verified_kintecus_too_slow_or_differs\BIOMD0000000316.kintpack -k "-ig:mass -INT:5 -obeymaxint"

THE COLUMNS OF SLOW_JOBS.csv
  order               run order (verified decks first, then by estimate)
  switches            recommended switches: the integrator that progressed fastest here
  estimate_here       run time extrapolated linearly from the furthest point reached here (Wine on a 2-core
                      cloud machine; a native Windows PC is usually faster). 'at least': the run was stopped
                      when this projection passed the time limit, so it is a lower bound. 'about': the run used
                      its whole time limit. 'unknown (stalls)': every integrator stopped producing output rows
                      for the whole stall window here; these may be very slow or truly stuck.
  estimate_basis      what the estimate is based on
  deck_verified       yes = the pack is confirmed correct; no / not checked = not confirmed
  has_sbml_reference  no = the run can finish but cannot be compared (no SBML reference)
  if_it_stalls_try    the other integrators to try (with --only ID --switches "...")
  previous_attempts   what each integrator did here

  Estimates here: 11 of 11 have one; 6 are estimated at under 1 hour and 5 more at
  1-24 hours. Linear extrapolation is rough: stiff systems often speed up or slow down part-way through.
  The Gear integrator (-INT:3) is deliberately not suggested.

IF A JOB STILL DOES NOT FINISH
  - No progress lines moving for a long time: stop it (Ctrl+C in the window) and rerun that pack with the next
    integrator in 'if_it_stalls_try'.
  - Very slow but moving: let it run; the progress line shows the remaining time it expects.
  - A FATAL error: send me the kintcheck_work\<ID>\ folder.
