BIOMD0000000554 (converter of 2026-11-07), switch species SW1-SW3 displayed.
Run:  Kintecus.exe -ig:mass -INT:1 -obeymaxint     -> CONC.TXT: SW2 stays 0 after t=200, SW3 stays 1 after 500,
                                                       SW1 stays 0 after 723; rows jump 197.5 -> 212.077 (no row at 200)
      Kintecus.exe -ig:mass -INT:5 -obeymaxint     -> SW2 = 1 from the 202.5 row, SW3 = 0 from 502.5: correct
                                                       (matches the SBML to 0.027%)
Small decks with the same three pair files and sampling 2.500000003 switch correctly on -INT:1
(slowdiag/jac/PF3, PF3c), so the deck size or stiffness matters.
