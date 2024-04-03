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

H1_13 = ham_ee(Sys1,[1 3]);
H2_13 = ham_ee(Sys2,[1 3]);
H1_12 = ham_ee(Sys1,[1 2]);
H2_12 = ham_ee(Sys2,[1 2]);
H1_23 = ham_ee(Sys1,[2 3]);
H2_23 = ham_ee(Sys2,[2 3]);

ok = areequal(H1_12,H2_12,1e-10,'abs') && ...
     areequal(H1_13,H2_13,1e-10,'abs') && ...
     areequal(H1_23,H2_23,1e-10,'abs');
