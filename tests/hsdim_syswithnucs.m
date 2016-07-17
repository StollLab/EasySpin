function [err,data] = test(opt,olddata)

% Test 2: spin system with nuclei
%================================================================
Sys = struct('S',1/2,'g',[2 3 4],'Nucs','14N','A',[6 34 4]);
v = hsdim(Sys);
err = ~areequal(v,6);
data = [];
