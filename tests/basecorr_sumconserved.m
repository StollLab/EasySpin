function [err,data] = test(opt,olddata)

%======================================================
% baseline + corrected data = data
%======================================================
N = 6345;
X = linspace(0,1,N);
Y = rand(1,N);
[c,b] = basecorr(Y,2,3);

err = any(abs(c+b-Y)>1e-10);
data = [];
