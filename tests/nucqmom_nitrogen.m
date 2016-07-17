function [err,data] = test(opt,olddata)

% 14N value
%======================================================
q = nucqmom('14N');
err = (q~=0.02044);
data = [];
