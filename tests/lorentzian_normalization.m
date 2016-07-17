function [err,data] = test(opt,olddata)

% Test whether Lorentzian line function is normalized
%======================================================
x0 = 34; fwhm = 2;
x = x0 + linspace(-100,100,1e5);
y = lorentzian(x,x0,fwhm);
Integral = trapz(x,y);

err = abs(Integral-1) > 1e-2;
data = [];
