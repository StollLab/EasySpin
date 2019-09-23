function [err,data] = test(opt,olddata)

% Lorentzian absorption and dispersion shapes
%======================================================
x = linspace(-1,1,100001);
x0 = 0;
fwhm = 0.2;

gamma = fwhm/sqrt(3); % distance between inflection point
pre = 2/pi/sqrt(3);

[yabs0,ydisp0] = lorentzian(x,x0,fwhm,0);
[yabs1,ydisp1] = lorentzian(x,x0,fwhm,1);
[yabs2,ydisp2] = lorentzian(x,x0,fwhm,2);

L = pre./(gamma+1i*sqrt(4/3)*(x-x0));

thr = 1e-3;
ok(1) = areequal(real(L),yabs0,thr,'abs') && areequal(-imag(L),ydisp0,thr,'abs');
ok(2) = areequal(deriv(x,yabs0),yabs1,thr,'abs') && areequal(deriv(x,ydisp0),ydisp1,thr,'abs');
ok(3) = areequal(deriv(x,yabs1),yabs2,thr,'abs') && areequal(deriv(x,ydisp1),ydisp2,thr,'abs');

err = any(~ok);
data = [];
