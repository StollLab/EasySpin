function [err,data] = test(opt,olddata)

% Test value
%======================================================
a = hbar;
b = planck/2/pi;
err = ~areequal(a,b,1e-12,'rel');
data = [];
