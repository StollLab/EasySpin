function ok = test(opt)

% Check that a dispersion signal with anisotropic widths can be computed with
% a Gaussian convolution width alone. The spectrum at mwPhase=pi/4 must be
% the phase-rotated mix of absorption and dispersion.

Sys.g = [2 2.05 2.1];
Sys.lwAnisoG = [50 100 150];
Sys.lw = [1 0];
Exp.mwFreq = 9.5;
Exp.Range = [315 345];
Exp.Harmonic = 0;

Exp.mwPhase = 0;
[x,spcAbs] = pepper(Sys,Exp);
Exp.mwPhase = pi/2;
[x,spcDisp] = pepper(Sys,Exp);
Exp.mwPhase = pi/4;
[x,spcMixed] = pepper(Sys,Exp);

ok = areequal(spcMixed,(spcAbs+spcDisp)/sqrt(2),1e-6,'rel') && ...
  max(abs(spcDisp))>0.1*max(spcAbs);

if opt.Display
  plot(x,spcAbs,x,spcDisp,x,spcMixed);
  legend('absorption','dispersion','\pi/4');
end

end
