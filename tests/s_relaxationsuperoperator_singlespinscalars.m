function ok = test()

% Scalar T1 and T2 for a spin-1/2: population difference relaxes with T1,
% coherences decay with T2, zero means no relaxation

H = diag([1 -1]);  % MHz, energy order opposite to product order
tol = 1e-12;

Sys.T1 = 1.5;
Sys.T2 = 0.5;

w = 1/(2*Sys.T1);  % rate per direction
RefT1 = zeros(4);
RefT1([1 4],[1 4]) = [w -w; -w w];
RefT2 = zeros(4);
RefT2(2,2) = 1/Sys.T2;
RefT2(3,3) = 1/Sys.T2;

Gamma = runprivate('s_relaxationsuperoperator',Sys,H);
ok(1) = areequal(Gamma,RefT1+RefT2,tol,'abs');

% Population difference decays with rate 1/T1
idx = [1 4];
k = eig(Gamma(idx,idx));
ok(2) = areequal(sort(k),[0; 1/Sys.T1],tol,'abs');

% No T1
Sys.T1 = 0;
Gamma = runprivate('s_relaxationsuperoperator',Sys,H);
ok(3) = areequal(Gamma,RefT2,tol,'abs');

% No T2
Sys.T1 = 1.5;
Sys.T2 = 0;
Gamma = runprivate('s_relaxationsuperoperator',Sys,H);
ok(4) = areequal(Gamma,RefT1,tol,'abs');
