function ok = test()

% 3 electron spins, check full vs simple

Sys.S = [1/2 1/2 1/2];
ee1 = [1 1 1];
ee2 = [2 2 2];
ee3 = [5 5 5];

Sys1.S = [1/2 1/2 1/2];
Sys1.ee = [ee1; ee2; ee3];
Sys2.S = Sys1.S;
Sys2.ee = [diag(ee1);diag(ee2);diag(ee3)];

H1 = ham_ee(Sys1);
H2 = ham_ee(Sys2);

ok = areequal(H1,H2,1e-10,'abs');
