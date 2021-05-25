function ok = test(opt)

% Simulations using Ci vs C2h

% Define spin system with C2h symmetry
Sys.g = [1.9,2,2.3];
Sys.Nucs = '1H';
Sys.A = [20 80 300];
Sys.AFrame = [-pi/4 -pi/2 0];
Sys.lw = 1;

SysPg = symm(Sys);
if ~strcmp(SysPg,'C2h')
  error('Incorrect symmetry %s of test system, should be C2h.',SysPg);
end

Exp.Range = [285 370];
Exp.mwFreq = 9.5;
Exp.Harmonic = 0;

Opt.GridSize = [50 3];
Opt.Method = 'perturb';

Opt.GridSymmetry = 'Ci';
[B,spc1] = pepper(Sys,Exp,Opt);
Opt.GridSymmetry = 'C2h';
[B,spc2] = pepper(Sys,Exp,Opt);

ok = areequal(spc1,spc2,2e-2,'rel');

if opt.Display
  title('Ci and C2h grids for C2h system');
  plot(B,spc1,B,spc2);
  legend('Ci','C2h');
end
