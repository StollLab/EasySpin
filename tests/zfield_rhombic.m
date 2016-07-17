function [err,data] = test(opt,olddata)

% Test of zfield() for a simple case
%================================================
Sys = struct('S',1,'D',[3 1]);
H0a = [1 0 1; 0 -2 0; 1 0 1];
H0b = zfield(Sys);

err = ~areequal(H0a,H0b);

data = [];
