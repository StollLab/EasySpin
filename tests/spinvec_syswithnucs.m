function [err,data] = test(opt,olddata)

Sys = struct('S',1/2,'g',[2 3 4],'Nucs','14N','A',[6 34 4]);
v = spinvec(Sys);
err = any(v~=[1/2 1]);
data = [];
