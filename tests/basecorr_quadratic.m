function [err,data] = test(opt,olddata)

%======================================================
% Quadratic baseline correction of quadratic data
%======================================================
N = 6345;
X = linspace(0,1,N);
Y = 1 - 0.234*X + 0.245*X.^2;
[c,b] = basecorr(Y,2,2);

err = any(abs(c)>1e-10);
data = [];
