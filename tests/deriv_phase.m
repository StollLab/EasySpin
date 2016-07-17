function [err,data] = test(opt,olddata)
% Test 2: correct phase
%================================================================
y = 1:100;
a = deriv(y);
err = any(a<=0);
data=[];
