function [err,data] = test(opt,olddata)

% Spin-orbit coupling between spins and orbital angular momenta can also be
% defined via electron-electron couplings between pure spins. The Hamiltonians
% constructed for these two cases must be identical.

rng_(3);
nSpins = randi_(3);
S = randi_(3,1,nSpins)/2;
L = randi_(2,1,nSpins);
soc = rand(nSpins,2);
orf = rand(nSpins,1);

% System containing S and L
SysSL.S = S;
SysSL.L = L;
SysSL.soc = soc;
SysSL.orf = orf;
if nSpins>1
  SysSL.ee = zeros(nchoosek(nSpins,2),1);
end

% System containing only spins
SysPureSpin.S = [S L];
len = 2*nSpins;
k = nchoosek(1:len,2);
eelen = nchoosek(len,2);
SysPureSpin.ee = zeros(eelen,1);
SysPureSpin.ee2 = zeros(eelen,1);
for m = 1:nSpins
  idx = (k(:,1)==m) & (k(:,2)==m+nSpins);
  SysPureSpin.ee(idx) = orf(m)*soc(m,1);
  SysPureSpin.ee2(idx) = orf(m)^2*soc(m,2);
end

% Calculate Hamiltonians
H1 = eeint(SysPureSpin);
H2 = soint(SysSL);

err = ~areequal(H1,H2,1e-10,'abs');
data = [];
