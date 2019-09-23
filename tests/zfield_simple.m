function [err,data] = test(opt,olddata)

% Test of zfield() for a simple case
%================================================
Sys = struct('S',3/2,'D',[-1 -1 2]);
H0a = diag([3 -3 -3 3]);
H0b = zfield(Sys);

err = ~areequal(H0a,H0b,1e-12,'abs');

data = [];
