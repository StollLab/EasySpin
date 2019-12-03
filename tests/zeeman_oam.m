function [err,data] = test(opt,olddata)

% Orbital angular momenta can also be defined as spins, therefore the two
% Hamiltonians should be identical.

rng_(5,'twister');

n = randi_(3);
S = randi_(3,1,n)/2;
L = randi_(3,1,n);

% System with S and L
SysSL.S = S;
SysSL.L = L;
SysSL.soc = zeros(n,1);
SysSL.g = rand(3*n,3);
SysSL.orf = rand(n, 1);
if n>1
  SysSL.ee = zeros(nchoosek(length(SysSL.S),2),1);
end

% System with [S L]
PureSpin.S = [SysSL.S,SysSL.L];
PureSpin.ee = zeros(nchoosek(length([SysSL.S,SysSL.L]),2),1);
PureSpin.g = [SysSL.g; zeros(3*n,3)];
for k = 1:n
  PureSpin.g(3*(n+k-1)+1:3*(n+k),:) = -diag(SysSL.orf(k)*ones(1,3));
end

% Compare Hamiltonians
H1 = zfield(PureSpin);
H2 = crystalfield(SysSL);

err = ~areequal(H1,H2,1e-10,'abs');
data = [];