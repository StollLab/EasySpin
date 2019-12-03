function err = test(opt,olddata)

% Test whether Lorentzian line function is normalized
%======================================================
x0 = 340;
fwhm = 20;
f = @(x) lorentzian(x,x0,fwhm);

n = 10;
xmin = x0 - n*fwhm;
xmax = x0 + n*fwhm;
Area = integral(f,xmin,xmax);

err = abs(Area-1) > 4e-2;
