function ok = test()


Sys = struct('S',1/2,'g',[2 2 2.2]);
Sys = nucspinadd(Sys,'14N',[3 3 9],[],[-1 -1 2]);
Sys = nucspinadd(Sys,'63Cu',[50 50 520]);
Sys = nucspinadd(Sys,'14N',[3 3 9],[],[-1 -1 2]);
A = [10 0 0; 0 15 0; 0 0 30];
Sys = nucspinadd(Sys,'14N',A);
Q = [1 0 0; 0 2 0; 0 0 3];
Sys = nucspinadd(Sys,'14N',[1],[],Q);

ok = numel(Sys.A)==45 & numel(Sys.Q)==45;