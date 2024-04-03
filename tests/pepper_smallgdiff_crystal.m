function ok = test(opt)

% Make sure simulations are accurate enough that differences between spectra
% differing minimally in their g values are accurate.

lwGL = [1 1];
SysA.g = 2.0001;
SysB.g = 2.0000;
SysA.lw = lwGL;
SysB.lw = lwGL;

Exp.mwFreq = 9.5;
Exp.Range = [320 360];
Exp.Harmonic = 0;
Exp.SampleFrame = [0 0 0];

Exp.nPoints = 100;
[B1,spcAcoarse] = pepper(SysA,Exp);
[B1,spcBcoarse] = pepper(SysB,Exp);
spcdiff1 = spcAcoarse-spcBcoarse;

Exp.nPoints = 20000;
[B2,spcAfine] = pepper(SysA,Exp);
[B2,spcBfine] = pepper(SysB,Exp);
spcdiff2 = spcAfine-spcBfine;

if opt.Display
  plot(B1,spcdiff1,B2,spcdiff2);
  legend('coarse','fine');
end

ok = areequal(max(spcdiff1),max(spcdiff2),0.15,'rel');
