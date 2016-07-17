function [err,data] = test(opt,olddata)

%======================================================
% Linear baseline correction of linear data
%======================================================

N = 1023;
Y = 1 + 0.234*linspace(0,1,N);
[c,b] = basecorr(Y,2,1);

err = any(abs(c)>1e-10);
data = [];
