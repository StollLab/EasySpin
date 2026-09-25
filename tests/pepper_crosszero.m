function ok = test()

% Spectrum over a symmetric field range across zero is even (absorption)
% or odd (first harmonic)

Sys.g = [2 2.1 2.2];
Sys.lw = 1;
Exp.mwFreq = 9.5;
Exp.Range = [-400 400];
Exp.nPoints = 2001;

Exp.Harmonic = 0;
[~,y0] = pepper(Sys,Exp);
Exp.Harmonic = 1;
[~,y1] = pepper(Sys,Exp);

% Each spectral point contains the integral up to the next point, so the
% mirrored spectrum is offset by one point
y0r = fliplr(y0);
y1r = -fliplr(y1);
ok = areequal(y0(1:end-1),y0r(2:end),1e-8,'rel') && ...
  areequal(y1(1:end-1),y1r(2:end),1e-8,'rel');
