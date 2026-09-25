function ok = test()

% Spectrum over an asymmetric field range across zero agrees with
% separate simulations over the positive part and the mirrored negative part

Sys.g = [2 2.1 2.2];
Sys.lw = 0.3;
Exp.mwFreq = 1;  % resonances between 32 and 36 mT
Exp.Harmonic = 0;

Exp.Range = [-40 100];
Exp.nPoints = 1401;
[~,y] = pepper(Sys,Exp);

Exp.Range = [0 100];
Exp.nPoints = 1001;
[~,ypos] = pepper(Sys,Exp);

Exp.Range = [0 40];
Exp.nPoints = 401;
[~,yneg] = pepper(Sys,Exp);

% Each spectral point contains the integral up to the next point, so the
% mirrored spectrum is offset by one point
yref = [fliplr(yneg(1:end-1)) ypos];

ok = areequal(y,yref,1e-3*max(y),'abs') && max(y(1:400))>0.1*max(y);
