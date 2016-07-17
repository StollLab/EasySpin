function [err,data] = test(opt,olddata)

% Test 2: Zero derivative at center
%======================================================
a = lorentzian(0,0,1.234,1);
err = abs(a)>1e-15;
data = [];
