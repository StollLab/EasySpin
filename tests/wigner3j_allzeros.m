function [err,data] = test(opt,olddata)

% Test 2
%======================================================
a = wigner3j(0,0,0,0,0,0);
b = 1;
err = ~areequal(a,b,0);
data = [];
