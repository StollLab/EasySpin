function ok = test(opt)

% Opt.ImmediateBinning (binning inside resfields_perturb) must give the same
% spectrum as the default interpolation/projection path.

Sys.g = [2.05 2.1 2.25];
Sys.Nucs = '63Cu';
Sys.A = [50 60 500];
Sys.lwpp = 2;

Exp.mwFreq = 9.5;
Exp.Range = [270 350];

Opt.Method = 'perturb';
Opt.GridSize = 91;

[B,spc0] = pepper(Sys,Exp,Opt);
Opt.ImmediateBinning = 1;
[B,spc1] = pepper(Sys,Exp,Opt);

if opt.Display
  plot(B,spc0,B,spc1);
  legend('default','ImmediateBinning');
end

ok = areequal(spc1,spc0,0.01,'rel');
