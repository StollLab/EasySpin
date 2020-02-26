function ok = test()

Sys = struct('S',[1/2 1/2],'g',[2 2 2;2.1 2.1 2.1],'ee',[10 10 10]);

eeint(Sys);
eeint(Sys,[1 2]);
H = eeint(Sys);
H = eeint(Sys,[1 2]);

ok = true;
