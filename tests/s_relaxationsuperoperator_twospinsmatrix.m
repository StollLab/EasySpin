function ok = test()

% Matrix T1 for several states: rate equation for the populations

tol = 1e-12;

% Spin-1 (3 states), relaxation only between states 1 and 2: the population
% difference of 1 and 2 relaxes with T1(1,2), the sum of their populations
% and the population of state 3 stay constant
H = diag([1 0 -1]);  % MHz
Sys.Spins = 1;
Sys.nElectrons = 1;
Sys.T1 = [0 2 0; 0 0 0; 0 0 0];
Sys.T2 = 0;
Gamma = runprivate('s_relaxationsuperoperator',Sys,H);
idxPop = [1 5 9];
x0 = [-1; 0; 1];  % deviation from equilibrium
t = 1.3;
x = expm(-Gamma(idxPop,idxPop)*t)*x0;
d = (x0(1)-x0(2))*exp(-t/Sys.T1(1,2));
s = x0(1)+x0(2);
ok(1) = areequal(x,[(s+d)/2; (s-d)/2; x0(3)],tol,'abs');

% Partially given lower triangle: missing entries are taken from the upper triangle
H = diag([1 2 3 4]);
N = 4;
T1 = [0 1 2 3; 0 0 4 5; 0 0 0 6; 0 0 0 0];
T1(3,1) = 7;  % only one reverse rate given
Sys.Spins = [1/2 1/2];
Sys.nElectrons = 2;
Sys.T1 = T1;
Sys.T2 = 0;
Gamma = runprivate('s_relaxationsuperoperator',Sys,H);
T1full = T1 + triu(T1,1).';
T1full(3,1) = 7;
W = 1./(2*T1full);
W(1:N+1:end) = 0;
K = W - diag(sum(W,1));
idxPop = 1:N+1:N^2;
ok(2) = areequal(Gamma(idxPop,idxPop),-K,tol,'abs');

% Diagonal elements of T1 are ignored
Sys.T1 = T1 + diag([1 1 1 1]);
GammaD = runprivate('s_relaxationsuperoperator',Sys,H);
ok(3) = areequal(GammaD,Gamma,tol,'abs');
