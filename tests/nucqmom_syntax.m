function [err,data] = test(opt,olddata)

% Test 1: syntax
%======================================================
q = nucqmom('14N');
q = nucqmom('63Cu,65Cu');
err = 0;
data = [];
