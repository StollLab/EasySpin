function ok = test()

% full hyperfine matrices, syntax check
rng(972);

% 1 electron, 2 nuclei
Sys.S = 1/2;
Sys.Nucs = '1H,1H';
Sys.A = [rand(3,3); rand(3,3)];
H1 = ham_hf(Sys);

% 2 electron, 1 nucleus
Sys.S = [1/2 1/2];
Sys.Nucs = '1H';
Sys.A = [rand(3,3), rand(3,3)];
Sys.ee = 1;
H2 = ham_hf(Sys);

% 2 electrons, 2 nuclei
Sys.S = [1/2 1/2];
Sys.Nucs = '1H,1H';
Sys.A = [rand(3,3), rand(3,3); rand(3,3) rand(3,3)];
Sys.ee = 1;
H3 = ham_hf(Sys);

% 1 electron, 1 spin-0 nucleus
Sys.S = 1/2;
Sys.Nucs = '12C';
Sys.A = rand(3,3);
H4 = ham_hf(Sys);

ok = true;
