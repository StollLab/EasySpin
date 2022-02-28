function ok = test(opt)

% Test Sys.Dstrain against an explicit loop over a distribution
% of D values.

Sys.S = 1;
Sys.g = 2;
Sys.D = 200;
Sys.DStrain = 50;
Sys.lwpp = 0.2;

Exp.mwFreq = 9.517;  % GHz
Exp.CenterSweep = [340 30];  % mT
Exp.Harmonic = 0;

Opt.GridSize = 121;
Opt.Method = 'perturb';

% Simulation using Sys.DStrain
[B,spc_DStrain] = pepper(Sys,Exp,Opt);

% Generate explicit distribution over D
D = Sys.D + linspace(-1,1,21)*1.5*Sys.DStrain;
w = gaussian(D,Sys.D,Sys.DStrain);
w = w/sum(w);

% Loop over D distribution, simulate spectra and add them
spc_explicit = 0;
Sys.DStrain = 0;
for k = 1:numel(D)
  Sys.D = D(k);
  spc_explicit = spc_explicit + w(k)*pepper(Sys,Exp,Opt);
end

ok = areequal(spc_DStrain,spc_explicit,0.03,'rel');

if opt.Display
  plot(B,spc_DStrain,B,spc_explicit)
end
