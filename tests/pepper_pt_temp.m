function ok = test(opt)

% temperature effects perturbation <-> matrix diagonalization

Sys.S = 1;
Sys.g = 2;
Sys.D = 400;
Sys.lwpp = 0.5;

Exp.mwFreq = 9.5;
Exp.CenterSweep = [330 60];
Exp.Temperature = 0.5;
Exp.Harmonic = 0;

Opt.Method = 'matrix';
[B,spc0] = pepper(Sys,Exp,Opt);
Opt.Method = 'perturb';
[~,spc1] = pepper(Sys,Exp,Opt);

if opt.Display
  plot(B,spc0,B,spc1);
  legend('matrix','perturb');
  legend boxoff
end

ok = areequal(spc0,spc1,2e-2,'rel');
