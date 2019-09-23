function [err,data] = test(opt,olddata)

% Make sure pepper returns orientations separately with Opt.Output = 'separate'
% in a crystal simulation

Sys.g = [2 2.1 2.2];
Sys.lwpp = 0.5;

ori1 = [0 0 0];
ori2 = [10 45 30]*pi/180;

Exp.mwFreq = 9.5;
Exp.Range = [300 400];

Opt.Output = 'separate';
Exp.CrystalOrientation = [ori1; ori2];
[B,spc3] = pepper(Sys,Exp,Opt);

Opt.Output = 'summed';
Exp.CrystalOrientation = ori1;
spc1 = pepper(Sys,Exp,Opt);
Exp.CrystalOrientation = ori2;
spc2 = pepper(Sys,Exp,Opt);

ok(1) = size(spc3,1)==2;
ok(2) = areequal([spc1;spc2],spc3,1e-10,'abs');

if opt.Display
  plot(B,spc3-[spc1;spc2]);
  xlabel('magnetic field (mT)');
end

err = ~ok;
data = [];