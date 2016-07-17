function [err,data] = test(opt,olddata)

%======================================================
% Linear baseline correction of linear data
%======================================================

N = 1023;
data = rand(N,N);

[c,baseline] = basecorr(data,2,0);
mean2 = mean(data,2);

err = any(abs(mean2-baseline(:,1))>1e-10);
data = [];
