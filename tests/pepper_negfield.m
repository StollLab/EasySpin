function ok = test()

% Spectrum at negative fields is the mirror image of the spectrum at
% positive fields

Sys.g = [2 2.1 2.2];
Sys.lw = 1;
Exp.mwFreq = 9.5;

Exp.Range = [290 350];
Exp.Harmonic = 0;
[~,y0p] = pepper(Sys,Exp);
Exp.Harmonic = 1;
[~,y1p] = pepper(Sys,Exp);

Exp.Range = [-350 -290];
Exp.Harmonic = 0;
[~,y0n] = pepper(Sys,Exp);
Exp.Harmonic = 1;
[~,y1n] = pepper(Sys,Exp);

% Each spectral point contains the integral up to the next point, so the
% mirrored spectrum is offset by one point
y0r = fliplr(y0p);
y1r = -fliplr(y1p);
ok = areequal(y0n(1:end-1),y0r(2:end),1e-8,'rel') && ...
  areequal(y1n(1:end-1),y1r(2:end),1e-8,'rel');
