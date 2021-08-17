function ok = test()

% Orbital angular momenta can also be defined as spins, therefore the two
% Hamiltonians should be identical.

rng(5);

n = randi(3);
S = randi(3,1,n)/2;
L = randi(3,1,n);
gL = rand(n, 1);
g = rand(3*n,3);

% System with S and L
SysSL.S = S;
SysSL.L = L;
SysSL.soc = zeros(n,1);
SysSL.g = g;
SysSL.orf = gL;
if n>1
  SysSL.ee = zeros(nchoosek(length(SysSL.S),2),1);
end

% System with [S L]
PureSpin.S = [SysSL.S,SysSL.L];
PureSpin.ee = zeros(nchoosek(length([SysSL.S,SysSL.L]),2),1);
PureSpin.g = [g; zeros(3*n,3)];
for k = 1:n
  gL = diag(SysSL.orf(k)*ones(1,3));
  PureSpin.g(3*(n+k-1)+1:3*(n+k),:) = gL;
end

% Compare Hamiltonians
H1 = zfield(PureSpin);
H2 = crystalfield(SysSL);

ok = areequal(H1,H2,1e-10,'abs');
