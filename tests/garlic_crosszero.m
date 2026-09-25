function ok = test()

% Spectrum over a symmetric field range across zero is odd (first harmonic)

Sys.g = 2;
Sys.Nucs = '14N';
Sys.A = 40;
Sys.lw = 1;
Exp.mwFreq = 9.5;
Exp.Range = [-400 400];
Exp.nPoints = 8001;

[~,y] = garlic(Sys,Exp);

ok = areequal(y,-fliplr(y),1e-10,'rel') && max(y)>0;
