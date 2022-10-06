function ok = test()

% Spin-orbit coupling between spins and orbital angular momenta can also be
% defined via electron-electron couplings between pure spins. The Hamiltonians
% constructed for these two cases must be identical.

rng(3);
nSpins = randi(3);
S = randi(3,1,nSpins)/2;
L = randi(2,1,nSpins);
soc = rand(nSpins,2);
gL = rand(nSpins,1);

% System containing S and L
SysSL.S = S;
SysSL.L = L;
SysSL.soc = soc;
SysSL.gL = gL;
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
  SysPureSpin.ee(idx) = soc(m,1);
  SysPureSpin.ee2(idx) = soc(m,2);
end

% Calculate Hamiltonians
H1 = ham_ee(SysPureSpin);
H2 = ham_so(SysSL);

ok = areequal(H1,H2,1e-10,'abs');
