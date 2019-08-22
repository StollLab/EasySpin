function [err,data] = test(opt,olddata)

% Test value
%======================================================
a = hartree;
b = 4.3597447222071e-18;
err = abs(a-b)/b > 1e-13;
data = [];
