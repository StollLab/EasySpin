% Slow-motion spectrum of copper phthalocyanine in sulfuric acid
%===============================================================================
clear, clf, clc

% CuPc parameters
%-------------------------------------------------------------------------------
CuPc.g = [2.0525 2.1994];
CuPc.Nucs = '63Cu,14N';
CuPc.n = [1 4]; % one 63Cu and four 14N
CuPc.A = [-54 -54 -608; 52.4 41.2 41.8]; % MHz

CuPc.tcorr = 10^-7.5; % seconds
CuPc.lw = 0.3; % mT

% Experimental parameters
%-------------------------------------------------------------------------------
Exp.mwFreq = 9.878; % GHz
Exp.CenterSweep = [325 100]; % mT

% Simulation options
%-------------------------------------------------------------------------------
% The following is the setting that tells chili to treat the 14N nuclei
% perturbationally using post-convolution, instead of including them into the
% full spin Hamiltonian matrix.
Opt.PostConvNucs = 2; % 2 is the index of 14N in the CuPc.Nucs
% Large basis set, since motion is slow, some settings are zero because of the
% axial symmetry of the electron+63Cu system
Opt.LLKM = [16 0 0 4];

% Simulation and plotting
%-------------------------------------------------------------------------------
chili(CuPc,Exp,Opt);
