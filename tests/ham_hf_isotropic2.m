function ok = test()

% expansion of isotropic hfc
% system with 2 electrons and 3 nuclei

aiso = [100 121; 34 56; 2 3];

Sys1.S = [1/2 1/2];
Sys1.ee = 1;
Sys1.Nucs = '1H,1H,1H';
Sys1.A = aiso;

Sys2.S = Sys1.S;
Sys2.ee = Sys1.ee;
Sys2.Nucs = Sys1.Nucs;
Sys2.A = kron(aiso,[1 1 1]);

H1 = ham_hf(Sys1);
H2 = ham_hf(Sys2);

ok = areequal(H1,H2,1e-10,'rel');
