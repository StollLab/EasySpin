function [err,data] = test(opt,olddata)

% Test transition selection for resfreqs_matrix

clear Sys Exp
Sys.S = 1;
Sys.D = -3*30e3*[1 0.000001];

Exp.Temperature = 0;

[nu, int] = resfreqs_matrix(Sys,Exp);

err = numel(nu)~=1;
data = [];
