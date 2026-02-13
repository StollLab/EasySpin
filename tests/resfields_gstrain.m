function ok = test()

% Test whether g-strain-based broadening gives correct line width

Sys.g = 2;
Sys.gStrain = 0.01;

Exp.mwFreq = 9.5;
Exp.Range = [300 400];

[B,I,W] = resfields(Sys,Exp);

ok = areequal(W/B,Sys.gStrain/Sys.g,1e-8,'rel');

end

