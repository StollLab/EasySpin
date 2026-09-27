function ok = test()

%=======================================================
% garlic should give the same integral intensity as pepper,
% for field sweeps and frequency sweeps
%=======================================================
Sys.g = 2.1;
Sys.Nucs = '1H';
Sys.A = 100;

% Field sweep
Sys.lwpp = 1;  % mT
Exp.mwFreq = 9.5;
Exp.Range = [300 360];
Exp.Harmonic = 0;

[x,y1] = pepper(Sys,Exp);
[x,y2] = garlic(Sys,Exp);
dx = x(2)-x(1);

ok(1) = areequal(sum(y1)*dx,sum(y2)*dx,0.001,'abs');

% Frequency sweep
Sys.lwpp = 5;  % MHz
clear Exp
Exp.Field = 340;
Exp.mwRange = [9.8 10.2];
Exp.Harmonic = 0;

[x,y1] = pepper(Sys,Exp);
[x,y2] = garlic(Sys,Exp);

ok(2) = areequal(sum(y1),sum(y2),1e-3,'rel');
