function ok = test()
clear;

Sys.S = [1/2 1/2];
Sys.g = [2.00 1.99];

Sys.J = -1; 
Sys.lwpp = 0.02; % mT
Sys.Nucs = '1H';
Sys.A = [2 0]; % MHz
Sys.Pop = 1/sqrt(2)*[0 1 -1 0];
Sys.PopMode = 'eigenbasis';

Expfield.CrystalOrientation = [0 pi/2];
Expfield.Range = [342 345.5]; % mT
Expfield.mwFreq = 9.6; % GHz
Expfield.Harmonic = 0;

Expfreq.CrystalOrientation = [0 pi/2];
Expfreq.Field = 344; % mT
Expfreq.mwRange = 9.6 + [-0.05 0.05]; % GHz
Expfreq.Harmonic = 0;

Opt.Method = 'hybrid';

[B,y1] = pepper(Sys,Expfield);

[~,y2] = pepper(Sys,Expfield,Opt);
Sys.lwpp = mt2mhz(Sys.lwpp); % mT
[~,y3] = pepper(Sys,Expfreq);

[~,y4] = pepper(Sys,Expfreq,Opt);


ok(1) = areequal(y1,y2,0.02,'rel');
ok(2) = areequal(y3,y4,0.02,'rel');



 