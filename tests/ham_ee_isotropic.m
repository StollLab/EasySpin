function ok = test()

% expansion of isotropic ee

ee = [100 121 37];

Sys1.S = [1/2 1/2 1/2];
Sys1.ee = ee;

Sys2.S = Sys1.S;
Sys2.ee = ee(:)*[1 1 1];

H1 = ham_ee(Sys1);
H2 = ham_ee(Sys2);

ok = areequal(H1,H2,1e-10,'abs');
