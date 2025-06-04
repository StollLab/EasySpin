function ok = test()

% Assert that rescaledata lsq mode works even if the two arrays are
% of different length

x0 = 0.7;
fwhm = 0.5;
scale = 132;
x = linspace(0,2,301);
y = gaussian(x,x0,fwhm);
xref = linspace(0,2,401);
yref = scale*gaussian(xref,x0,fwhm);

[~,scale_fit] = rescaledata(y,yref,'lsq');

ok = areequal(scale_fit,scale,1e-5,'rel');
