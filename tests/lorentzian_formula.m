function [err,data] = test(opt,olddata)

% Test: Check against formula (no derivative)
%======================================================
x0 = 3; fwhm = 2;
x = -3:6;
y = lorentzian(x,x0,fwhm,0);
gam = fwhm/sqrt(3);
xi = (x-x0)/gam;
y0 = 2/pi/sqrt(3)/gam./(1+4/3*xi.^2);

err = any(abs(y-y0)/y>1e-10);
data = [];
