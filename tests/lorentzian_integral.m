function [err,data] = test(opt,olddata)

% Test 1: Positive values, integral normalized
%======================================================
x = 10000*linspace(-1,1,1e5); x0 = 34; w = 10;
y = lorentzian(x,x0,w);
Area = sum(y(1:end-1).*diff(x));
Area = trapz(x,y);

err = (abs(Area-1)>0.01) | any(y<=0);
data = [];
