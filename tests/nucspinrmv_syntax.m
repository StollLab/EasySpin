function ok = test()

Sys = struct('S',1/2,'g',[2 2.1 2.2]);
Sys = nucspinadd(Sys,'14N',[1 2 3]);
Sys = nucspinadd(Sys,'1H',[23 34 45]);
nucspinrmv(Sys,1);
nucspinrmv(Sys,2);
nucspinrmv(Sys,[1 2]);

ok = true;
