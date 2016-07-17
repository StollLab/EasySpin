function [err,data] = test(opt,olddata)

% Test 1: syntax
%======================================================
an = rand(1,3);
R1 = erot(an);
eulang(R1);
[a,b,c] = eulang(R1);
aa = eulang(R1);
err = 0;
data = [];
