function ok = test(opt)

% Test whether Gaussian, Lorentzian, and Gaussian+Lorentzian broadenings always
% give requested harmonic.

Sys.g = 1.9995;
Exp.mwFreq = 9.5;
Exp.Range = [339 340];
lwppG = 0.1;
lwppL = 0.08;

isfirstharmonic = @(spc) areequal(max(spc),-min(spc),1e-3,'rel');

% Gaussian only
Sys.lwpp = lwppG;
[B,spc1] = pepper(Sys,Exp);
ok(1) = isfirstharmonic(spc1);

% Gaussian with zero Lorentzian
Sys.lwpp = [lwppG 0];
[B,spc2] = pepper(Sys,Exp);
ok(2) = isfirstharmonic(spc2);

% Lorentzian with zero Gaussian
Sys.lwpp = [0 lwppL];
[B,spc3] = pepper(Sys,Exp);
ok(3) = isfirstharmonic(spc3);

% Combined Gaussian and Lorentzian
Sys.lwpp = [lwppG lwppL];
[B,spc4] = pepper(Sys,Exp);
ok(4) = isfirstharmonic(spc4);

if opt.Display
  subplot(2,2,1)
  plot(B,spc1)
  title('1: Gaussian only')
  subplot(2,2,2)
  plot(B,spc2)
  title('2: Gaussian with zero Lorentzian')
  subplot(2,2,3)
  plot(B,spc3)
  title('3: Lorentzian with zero Gaussian')
  subplot(2,2,4)
  plot(B,spc4)
  title('4: Gaussian and Lorentzian')
end
