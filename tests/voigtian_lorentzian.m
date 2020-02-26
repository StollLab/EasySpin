function ok = test()

% Compare voigtian with Lorentzian width 0 to explicit Gaussian

x = linspace(0,1,101);

x0 = 0.4;
fwhm = [0 0.1];

y1 = voigtian(x,x0,fwhm);
y2 = lorentzian(x,x0,fwhm(2));

ok = areequal(y1,y2,1e-8,'abs');
