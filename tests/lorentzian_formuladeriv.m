function [err,data] = test(opt,olddata)

% Test: Check against formula (first derivative)
%======================================================
x0 = 3;
x = -3:6;
fwhm = 2;
gam = fwhm/sqrt(3);

y = lorentzian(x,x0,fwhm,1);
y0 = -(16/3/sqrt(3)/pi)*1/gam^2*((x-x0)/gam).*(1+4/3*((x-x0)/gam).^2).^-2;

err = any(abs(y-y0)/y>1e-10);
data = [];
