function [err,data] = test(opt,olddata)

%======================================================
% Test 3: several fields and several nuclei
%======================================================
B = 100:150;
nu = larmorfrq('14N,15N,63Cu',100:150);
err = any(size(nu)~=[numel(B) 3]);
data = [];
