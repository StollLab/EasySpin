function ok = test(opt)

% The dispersion spectrum of a frequency sweep is the mirror image of the
% one of the field sweep: a spin with resonance field B in a field sweep at
% mwFreq0 resonates at frequency mwFreq0*Field0/B in a frequency sweep at Field0.

Sys.g = [2.01 2.00];
Sys.tcorr = 1e-9;

mwFreq0 = 9.8;  % GHz
Field0 = 349;  % mT

ExpB.mwFreq = mwFreq0;
ExpB.Range = [345 353];
ExpB.Harmonic = 0;
ExpB.mwPhase = pi/2;
[B,yB] = chili(Sys,ExpB);

ExpF.Field = Field0;
ExpF.mwRange = mwFreq0*Field0./fliplr(ExpB.Range);
ExpF.Harmonic = 0;
ExpF.mwPhase = pi/2;
[nu,yF] = chili(Sys,ExpF);

% Field sweep spectrum mapped onto the frequency axis
yBF = interp1(B,yB,mwFreq0*Field0./nu,'linear',0);

c = corrcoef(yF,yBF);
ok = c(1,2)>0.9;

if opt.Display
  plot(nu,yF/max(abs(yF)),nu,yBF/max(abs(yBF)));
  legend('frequency sweep','field sweep, mapped');
  title(sprintf('correlation %g',c(1,2)));
end

end
