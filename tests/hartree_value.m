function [err,data] = test(opt,olddata)

% Test value
%======================================================
a = hartree;
b = 4.359744650e-18;
err = abs(a-b)/b > 1e-12;
data = [];
