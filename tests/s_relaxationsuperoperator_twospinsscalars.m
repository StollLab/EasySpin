function ok = test()

% Scalar T1 and T2 for two spins-1/2 (4 states): same as a matrix with
% equal values, trace is conserved, coherences decay with T2 (two electrons:
% all coherences are electron coherences; electron and nucleus: nuclear
% coherences do not decay)

H = diag([1 2 3 4]);  % MHz
tol = 1e-12;
N = 4;
idxPop = 1:N+1:N^2;
idxCoh = setdiff(1:N^2,idxPop);

Sys.Spins = [1/2 1/2];
Sys.nElectrons = 2;
Sys.T1 = 1.5;
Sys.T2 = 0.5;
Gamma = runprivate('s_relaxationsuperoperator',Sys,H);

% Same as matrix with equal values
SysM = Sys;
SysM.T1 = triu(1.5*ones(N),1);
SysM.T2 = triu(0.5*ones(N),1);
GammaM = runprivate('s_relaxationsuperoperator',SysM,H);
ok(1) = areequal(Gamma,GammaM,tol,'abs');

% Populations: trace conserved, all other modes decay with rate N/(2*T1)
GammaPop = Gamma(idxPop,idxPop);
ok(2) = areequal(sum(GammaPop,1),zeros(1,N),tol,'abs');
k = sort(eig(GammaPop));
ok(3) = areequal(k,[0; N/(2*Sys.T1)*ones(N-1,1)],tol,'abs');

% Coherences: decay with T2, no coupling to populations
ok(4) = areequal(Gamma(idxCoh,idxCoh),eye(N^2-N)/Sys.T2,tol,'abs');
ok(5) = ~any(Gamma(idxPop,idxCoh),'all') && ~any(Gamma(idxCoh,idxPop),'all');

% Electron and nucleus: nuclear coherences (product states 1-2 and 3-4) do
% not decay, electron coherences decay with T2
Sys.nElectrons = 1;
Gamma = runprivate('s_relaxationsuperoperator',Sys,H);
R2 = [0 0 1 1; 0 0 1 1; 1 1 0 0; 1 1 0 0]/Sys.T2;
ok(6) = areequal(diag(Gamma(idxCoh,idxCoh)),R2(idxCoh).',tol,'abs');
