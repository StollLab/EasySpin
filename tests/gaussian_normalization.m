function err = test(opt,olddata)

% Test whether Gaussian linefunction is normalized
%======================================================
x0 = 340;
fwhm = 20;
f = @(x) gaussian(x,x0,fwhm);

n = 3;
xmin = x0-n*fwhm;
xmax = x0+n*fwhm;
Area = integral(f,xmin,xmax);

err = abs(Area-1) > 1e-10;
