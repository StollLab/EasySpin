function ok = test()

Sys = struct('S',1/2,'g',[2 2.1 2.2]);
Sys = nucspinadd(Sys,'14N',[1 2 3]);
Sys = nucspinadd(Sys,'1H',[23 34 45]);

Sys1  = nucspinkeep(Sys,1);
Sys2  = nucspinkeep(Sys,2);
Sys12 = nucspinkeep(Sys,[1 2]);

ok = size(Sys1.A,1)==1 && size(Sys2.A,1)==1 && size(Sys12.A,1)==2;
