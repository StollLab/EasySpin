function ok = test()

x = linspace(0,1,1002);
x0 = 0.4;
fwhm = [0.1 0.07];

y = voigtian(x,x0,fwhm);
y = voigtian(x,x0,fwhm,0);
y = voigtian(x,x0,fwhm,1);
y = voigtian(x,x0,fwhm,2);

ok = true;
