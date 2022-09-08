function ok = test()

% explicit comparison of hyperfine Hamiltonian

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A = [10 12 -44];
Sys.g = 2;

H0 = [-11 0 0 -0.5; 0 11 5.5 0; 0 5.5 11 0; -0.5 0 0 -11];
H1 = ham_hf(Sys);

ok = areequal(H1,H0,1e-10,'rel');
