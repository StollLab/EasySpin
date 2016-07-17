function [err,data] = test(opt,olddata)

%======================================================
% Test 2: several fields
%======================================================
Nuc = '14N';
B = 100:200;
nu1 = larmorfrq(Nuc,B);
nu2 = abs(nucgval(Nuc))*B*nmagn/planck/1e9;
err = any(abs(nu1-nu2.')>1e-10);
data = [];
