STATUS: FIXED in the 2026-11-23e build (md5 88f7fe5e...). Rerun 2026-11-24 at t=10 on -INT:1, -INT:2, -INT:4 and
-INT:5: X = 50, Y = 50, A = 0.832013, Z = 6.81989 (all exact; -INT:5 writes its last row at t=10.0001).

TIME in USER functions - build compiled Nov 12 2026 18:39:31

Run:  Kintecus.exe -ig:mass -INT:1 -obeymaxint      (same result with -INT:5)

Expected (t = run time in seconds):
  X : dX/dt = TIME                      -> X = t^2/2
  Y : dY/dt = [TIME_], the old way       -> Y = t^2/2     (source==>TIME_ at k=1 makes [TIME_] = t)
  A : dA/dt = -0.1*sin(TIME)*A           -> A = exp(-0.1*(1-cos t))
  Z : dZ/dt = A*TIME/(1+TIME)

Got: Y is correct (0.053749 at t=0.327869 = t^2/2), but X stays 1e-55, A stays 1.0 and Z stays 1e-55:
TIME evaluates to 0 in every line that lists it.

Other forms tried (single-line decks, -INT:1 and -INT:5), all with TIME = 0:
  S[USER;TIME:(TIME)]         X stays 1e-55
  S[USER;TIME:(1+TIME)]       X = t        (so TIME = 0)
  S[USER;cA,TIME:(cA*TIME)]   X stays 1e-55
  TIME also declared as a species row in SPECIES.DAT, with S[USER;TIME:(1+TIME)] or S[USER;cTIME:(1+cTIME)]:
                              X = t
  S[USER;cA:(TIME+0*cA)]      (TIME used but not listed) X = -t, i.e. TIME reads -1
