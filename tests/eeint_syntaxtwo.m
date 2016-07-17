function [err,data] = test(opt,olddata)

%======================================================
% Test 1: Syntax, 2 electrons
%======================================================
Sys = struct('S',[1/2 1/2],'g',[2 2 2;2.1 2.1 2.1],'ee',[10 10 10]);

eeint(Sys);
eeint(Sys,[1 2]);
H = eeint(Sys);
H = eeint(Sys,[1 2]);
err = 0;
data = [];
