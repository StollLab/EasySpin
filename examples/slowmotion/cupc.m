% Slow-motion spectrum of copper phthalocyanine in sulfuric acid
%==========================================================================
clear, clf

% CuPc parameters
%--------------------------------------------------------------------------
g = [2.0525 2.1994];
CuPc = struct('g',g,'Nucs','63Cu,14N','n',[1 4]);
CuPc.A = [-54 -54 -608; 52.4 41.2 41.8];

CuPc.logtcorr = -7.35;
CuPc.lw = 0.3;

% Experimental parameters
%----------------------------------------------------
Exp = struct('mwFreq',9.878,'CenterSweep',[325 100]);
Exp.nPoints = 1e4;

% Simulation and plotting
%----------------------------------------------------
Opt.LLKM = [36 30 0 8];
chili(CuPc,Exp,Opt);
