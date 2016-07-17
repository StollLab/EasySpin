function [err,data] = test(opt,olddata)

%======================================================
% Empty isotope list
%======================================================
nu1 = larmorfrq('',150);
nu2 = larmorfrq('1H',[]);
err = ~isempty(nu1) | ~isempty(nu2);
data = [];
