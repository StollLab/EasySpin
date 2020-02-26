function ok = test()

x = linspace(0,1,101);
x0 = 0.4;
fwhm = [0.1 0.07];

y = voigtian(x,x0,fwhm);
[ymax,idx] = max(y);

ok = (idx==41);
