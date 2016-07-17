% CF3 radical in gas matrix at low temperature
%==========================================================================
% see Edlund et al, Molecular Physics 32 (1976) 49
clear, clf

% Spin system with 3 19F nuclei
%--------------------------------------------------------------------
Sys.Nucs = '19F,19F,19F';
Sys.lwpp = 0.4;
Sys.g = [2.0042, 2.0024];

% Hyperfine coupling tensors
%--------------------------------------------------------------------
% 19F hyperfine coupling tensors: identical principal values
A = [80 87 262]*2.8; % in MHz
Sys.A = [A; A; A];
% ...but different orientations
beta = 17.8;
Sys.AFrame = [0 beta 0; -120 beta 0; 120 beta 0]*pi/180;

% Other parameters
%--------------------------------------------------------------------
Exp.mwFreq = 9.27; Exp.Range = [285 375];

% use perturbation theory instead of matrix diagonalization
% this is often useful for systems with several nuclei
Opt.Method = 'perturb2';

% Simulation and graphical rendering
%--------------------------------------------------------------------
pepper(Sys,Exp,Opt);

title(sprintf('CF3 radical, randomly oriented, %g GHz',Exp.mwFreq));
