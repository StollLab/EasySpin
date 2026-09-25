function ok = test()

% Spectrum at negative fields is the mirror image of the spectrum at
% positive fields (first harmonic, so the sign changes as well)

Sys.g = 2;
Sys.Nucs = '14N';
Sys.A = 40;
Sys.lw = 0.1;
Exp.mwFreq = 9.5;

Exp.Range = [330 350];
[~,yp] = garlic(Sys,Exp);
Exp.Range = [-350 -330];
[~,yn] = garlic(Sys,Exp);

ok = areequal(yn,-fliplr(yp),1e-10,'rel');
