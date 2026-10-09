function ok = test()

% Non-symmetric full A tensor: perturbation theory agrees with exact diagonalization
% of the electron Zeeman + hyperfine Hamiltonian (no nuclear Zeeman, as in perturbation theory)

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.g = 2;
Sys.A = [10 30 -20; -15 5 25; 40 -10 20];  % MHz

Exp.Field = 12000;  % mT, high field, field along z

Opt.PerturbOrder = 2;
Pos1 = sort(resfreqs_perturb(Sys,Exp,Opt));

H = ham_ez(Sys,[0 0 Exp.Field]) + ham_hf(Sys);
E = eig((H+H')/2);  % ascending; levels 1,2: mS=-1/2 (mI=+1/2,-1/2), levels 3,4: mS=+1/2 (mI=-1/2,+1/2)
Pos2 = sort([E(3)-E(2); E(4)-E(1)]);

ok = areequal(Pos1,Pos2,1e-8,'rel');
