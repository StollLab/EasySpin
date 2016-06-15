function [err,data] = test(opt,olddata)

%======================================================
% Test 1: 1 value
%======================================================
Nuc = '1H';
B = 100;
nu1 = larmorfrq(Nuc,B);
nu2 = abs(nucgval(Nuc))*B*nmagn/planck/1e9;
err = abs(nu1-nu2)>1e-10;
data = [];
