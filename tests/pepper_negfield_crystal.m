function ok = test()

% Crystal with several sites at negative fields: resonance fields are the
% mirror images of those at positive fields

Sys.g = [2 2.1 2.2];
Sys.gFrame = [0.3 0.4 0.5];
Exp.mwFreq = 9.5;
Exp.CrystalSymmetry = 'C2h';
Exp.SampleFrame = [0.1 0.2 0.3];

Exp.Range = [290 350];
Pp = resfields(Sys,Exp);
Exp.Range = [-350 -290];
Pn = resfields(Sys,Exp);

Exp.Range = [290 350];
[~,yp] = pepper(Sys,Exp);
Exp.Range = [-350 -290];
[~,yn] = pepper(Sys,Exp);

ok = areequal(Pn,-Pp,1e-10,'abs') && areequal(max(abs(yn)),max(abs(yp)),1e-6,'rel');
