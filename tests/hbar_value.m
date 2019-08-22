function [err,data] = test(opt,olddata)

% Test value
%======================================================
a = hbar;
b = planck/2/pi;
err = ~areequal(a,b);
data = [];
