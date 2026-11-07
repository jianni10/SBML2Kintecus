Kintecus V2027 (build 2026-10-30 22:51:25): a species driven below zero stalls the integration
=========================================================================================

Deck: NEG_FROM_ZERO_stall/ (MODEL.DAT, SPECIES.DAT, PARM.DAT)

    -0.8275  0 0  source==>A                         (negative rate constant = constant drain of A)
     0.678   0 0  B==>sink
     1       0 0  source+S[USER;cB:(0.662*cB)]==>A   (production of A from B)
    A(0) = 0, B(0) = 1.25, source and sink constant 1.0; 0 to 10 s.

At t = 0 the two terms cancel exactly (0.662 x 1.25 = 0.8275), then B decays and A goes negative:
    A(t) = -0.8275 t + 0.8275 (1 - exp(-0.678 t)) / 0.678          A(10) = -7.0559 (exact)

  Kintecus -ig:mass -INT:1 -obeymaxint -i      no output after t = 0 (stalls; same with -INT:2 and -INT:5)
  Kintecus -ig:mass -INT:5:10:0 -obeymaxint -i reaches t = 10 with A = -7.05589 (correct), but then loops on
                                               "----AFTER DASPK, IDID= -33 / ----TEMP, RTEMP= 9.99999999999999
                                               10.0000000000000 / There is a possible problem with your input
                                               to DASPK." and never writes KDONE.OUT
  The same deck with -0.80 instead of -0.8275 (A starts upward, then crosses zero) runs with -INT:1 (A(10) =
  -6.78) but stalls with -INT:5.

Why it matters for the BioModels packs: many SBML models have state variables that legitimately go to or
below zero (SBML parameters with rate rules, e.g. deviations from a basal level; in 12 of the 79 slow
models the libRoadRunner reference itself has variables going negative). Written exactly ("signed" lines,
as above), they need Kintecus to integrate through zero; -INT:5:10:0 (DASPK without the non-negativity
constraint) does, apart from the end-of-run loop.
