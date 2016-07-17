function [err,data] = test(opt,olddata)

%======================================================
% Test 2: Syntax, more than 2 electron spins
%======================================================
Sys = struct('S',[1/2 1/2 1/2],'g',[2 2 2; 2.1 2.1 2.1; 2.2 2.2 2.2]);
Sys.ee = [1 1 1; 2 2 2; 3 3 3];
H = eeint(Sys);
H = eeint(Sys,[1 2]);
H = eeint(Sys,[1 3]);
H = eeint(Sys,[2 3]);
err = 0;
data = [];
