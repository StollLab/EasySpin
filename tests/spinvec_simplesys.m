function [err,data] = test(opt,olddata)

% Test 1
%================================================================
Sys = struct('S',1,'g',[3 4 5]);
v = spinvec(Sys);
err = (v~=1);
data = [];
