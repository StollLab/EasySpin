function ok = test()

Sys = struct('S',1/2,'g',[2 2.1 2.2]);

S = nucspinadd(Sys,'14N',[1 2 3]);
S = nucspinadd(Sys,'14N',[1 2 3],[]);
S = nucspinadd(Sys,'14N',[1 2 3],[],[]);
S = nucspinadd(Sys,'14N',[1 2 3],[],[6 5 4],[1 2 3]);

ok = true;
