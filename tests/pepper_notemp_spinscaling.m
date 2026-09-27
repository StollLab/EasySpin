function ok = test()

% Without temperature, integrated intensities follow the high-temperature
% limit: proportional to sum of S(S+1) over all electrons

Exp.mwFreq = 9.5;
Exp.Range = [320 360];
Exp.Harmonic = 0;

Sys.g = 2;
Sys.lwpp = 1;

Sys.S = 1/2;
[~,y0] = pepper(Sys,Exp);

Sys.S = 1;
Sys.D = 100;
[~,y1] = pepper(Sys,Exp);

Sys = rmfield(Sys,'D');
Sys.S = [1/2 1/2];
Sys.g = [2; 2];
Sys.ee = 20;
[~,y2] = pepper(Sys,Exp);

ok(1) = areequal(sum(y1)/sum(y0),8/3,0.01,'rel');
ok(2) = areequal(sum(y2)/sum(y0),2,0.01,'rel');
