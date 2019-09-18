function [err,data] = test(opt,olddata)

%=======================================================
% Explicit field-sweep doublet, fast method
%=======================================================

Sys.S = 1/2;
Sys.g = [2.01 2.005 2.002];
Sys.tcorr = 10e-9;

Exp.mwFreq = 9.5;
Exp.Range = [336 341];
Exp.Harmonic = 0;

Opt.Rescale = true;
[x,y1] = chili(Sys,Exp,Opt);
Opt.Rescale = false;
[x,y2] = chili(Sys,Exp,Opt);

err = ~areequal(y1,y2,1e-10*max(y1));
data = [];

if opt.Display
  plot(x,y1,x,y2);
  legend('Rescale=true','Rescale=false');
end
