function [err,data] = test(opt,olddata)

%======================================================
% Test 1: Simple case 1 nucleus
%======================================================
Sys = struct('S',1/2,'g',2,'Nucs','2H','A',[1 2 3]);
Sys.Q = [4 6 8];
QQ([1 15 22 36]) = 13;
QQ([8 29]) = 10;
QQ([3 13 24 34]) = -1;
HQ = nquad(Sys);

err = any(abs(HQ(:)-QQ(:))>1e-10);
data = [];
