function ok = test(opt)

% Check that field modulation works without convolutional line width: the
% modulated spectrum is the modulation-broadened stick spectrum, and the
% default harmonic is 1.

Sys.g = 2;
Exp.mwFreq = 9.5;
Exp.Range = [335 345];

Exp.Harmonic = 0;
[x,spc0] = pepper(Sys,Exp);

Exp = rmfield(Exp,'Harmonic');
Exp.ModAmp = 0.5;
[x,spc1] = pepper(Sys,Exp);

spcRef = fieldmod(x,spc0,Exp.ModAmp,1);

ok = areequal(spc1,spcRef,1e-10,'rel');

if opt.Display
  plot(x,spcRef,x,spc1,'--');
  legend('fieldmod of stick spectrum','pepper');
end

end
