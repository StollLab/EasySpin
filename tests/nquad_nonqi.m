function [err,data] = test(opt,olddata)

%======================================================
% Test 2: no quadrupole interaction
%======================================================
Sys = struct('S',1/2,'g',[2 2 2],'Nucs','1H','A',[1 2 3]);
HQ = nquad(Sys);

err = any(HQ(:));
data = [];
