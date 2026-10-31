function ok = test()

% Relaxation times in matrices refer to product states, independent of the
% energy order of the eigenstates and of the basis returned by eig for
% degenerate eigenstates

tol = 1e-10;
N = 4;

% T2 only for the coherence between product states 1 and 3, T1 only between
% product states 2 and 4
Sys.T1 = zeros(N); Sys.T1(2,4) = 1;
Sys.T2 = zeros(N); Sys.T2(1,3) = 0.2;

% Reference: diagonal Hamiltonian with energies in product order
Ref = runprivate('s_relaxationsuperoperator',Sys,diag([4 3 2 1]));

% Energy order different from product order
Gamma = runprivate('s_relaxationsuperoperator',Sys,diag([2 4 1 3]));
ok(1) = areequal(Gamma,Ref,tol,'abs');

% Fully degenerate Hamiltonian with tiny perturbation: eig returns an
% arbitrary basis
rng(1);
X = rand(N) + 1i*rand(N);
[Gamma,minOverlap] = runprivate('s_relaxationsuperoperator',Sys,1e-13*(X+X'));
ok(2) = areequal(Gamma,Ref,tol,'abs');
ok(3) = abs(minOverlap-1)<tol;

% Weakly mixed states: overlap below 1, superoperator close to reference
H = diag([4 3 2 1]);
H(1,2) = 0.01; H(2,1) = 0.01;
[Gamma,minOverlap] = runprivate('s_relaxationsuperoperator',Sys,H);
ok(4) = minOverlap<1 && minOverlap>0.99;
ok(5) = areequal(Gamma,Ref,0.05,'rel');
