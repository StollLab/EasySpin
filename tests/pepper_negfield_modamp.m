function ok = test()

% Field modulation at negative fields: first-harmonic spectrum is the
% inverted mirror image of the spectrum at positive fields

Sys.g = [2 2.1 2.2];
Sys.lw = 1;
Exp.mwFreq = 9.5;
Exp.ModAmp = 0.5;

Exp.Range = [290 350];
[~,yp] = pepper(Sys,Exp);
Exp.Range = [-350 -290];
[~,yn] = pepper(Sys,Exp);

% Each spectral point contains the integral up to the next point, so the
% mirrored spectrum is offset by one point
yr = -fliplr(yp);
ok = areequal(yn(2:end-2),yr(3:end-1),1e-6,'rel');
