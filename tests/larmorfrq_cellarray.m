function [err,data] = test(opt,olddata)

%======================================================
% Test 6: cell array of isotopes
%======================================================
nu = larmorfrq({'1H','14N'},300);
err = (numel(nu)~=2);
data = [];
