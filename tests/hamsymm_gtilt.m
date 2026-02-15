function ok = test()

Sys.g = [2 2.1 2.2];
Sys.gFrame = deg2rad([10 20 30]);

[pg,R] = hamsymm(Sys);
symmFrame = eulang(R.');

ok = areequal(symmFrame,Sys.gFrame,1e-10,'rel');
