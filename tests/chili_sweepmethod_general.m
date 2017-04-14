function [err,data] = test(opt,olddata)

%================================================================
% Compare approximate and explicit field-swept spectra (general)
%================================================================

Sys.g = [2.01 2.003];
Sys.Nucs = '14N';
Sys.A = [20 20 100];
Sys.tcorr = 1e-9;

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.nPoints = 200;

Opt.LiouvMethod = 'general';

[x,y1] = chili(Sys,Exp,Opt);

Opt.ExplicitFieldSweep = true;
[x,y2] = chili(Sys,Exp,Opt);

if opt.Display
  plot(x,y1,x,y2);
  legend('approximate','explicit');
end

err = ~areequal(y1,y2,1e-3*max(y1));

data = [];