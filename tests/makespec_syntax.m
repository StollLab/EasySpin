function [err,data] = test(opt,olddata)

% Test 1: syntax
%====================================================
s1 = makespec([0 1],1000,0.5,1);
[x,s2] = makespec([0 1],1000,0.5,1);
err = any(abs(s1-s2))>1e-10;
data = [];
