function ok = test()

% Test D (2 variants) and B2 terms against explicit formula

S = 1;
Sys1.S = S;
Sys1.g = 2;
Sys2 = Sys1;
Sys3 = Sys1;

% Zero-field splitting parameters
D = 100*rand; E = D*rand;
Ddiag = [-1,-1,2]/3*D + [+1,-1,0]*E;

% (0) Explicit formula -------------------------------------------
[Sx,Sy,Sz] = sop(S,'x','y','z');
H0 = Ddiag(1)*Sx^2 + Ddiag(2)*Sy^2 + Ddiag(3)*Sz^2;

% (1) D and E ----------------------------------------------------
Sys1.D = [D,E];
H1 = ham_zf(Sys1);

% (2) D diagonal -------------------------------------------------
Sys2.D = Ddiag;
H2 = ham_zf(Sys2);

% (3) B2 parameters ----------------------------------------------
Sys3.B2 = [E 0 D/3 0 0];
H3 = ham_zf(Sys3);

% Compare all Hamiltonians
threshold = 1e-7;
ok(1) = areequal(H0,H1,threshold,'abs');
ok(2) = areequal(H0,H2,threshold,'abs');
ok(3) = areequal(H0,H3,threshold,'abs');
