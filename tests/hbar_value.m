function [err,data] = test(opt,olddata)

% Test value
%======================================================
a = hbar;
b = 1.054571800e-34;
err = abs(a-b)/b > 1e-12;
data = [];
