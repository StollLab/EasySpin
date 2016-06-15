function [err,data] = test(opt,olddata)

% Test 1: spin system
%================================================================
Sys = struct('S',1,'g',[3 4 5]);
v = hsdim(Sys);
err = (v~=3);
data = [];
