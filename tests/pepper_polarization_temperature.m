function ok = test()
clear

Sys.S = 1/2;
Sys.g = 2;
Sys.lwpp = .1;

Exp.nPoints = 1024;
Exp.mwFreq = 9.5;
Exp.Temperature = 0.1;

[~,y1] = pepper(Sys,Exp);

Sys.PopMode = 'highfield';
BoltzmannPreFactor = -1e6*planck/boltzm/Exp.Temperature;
Sys.Pop = exp(BoltzmannPreFactor*[0 9500]);

[~,y2] = pepper(Sys,Exp);

ok = areequal(y1,y2,1e-10,'rel');



