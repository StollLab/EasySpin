function ok = test()

fwhm = 1.1; % full width at half maximum

% Calculate Lorentzian lineshape
x0 = 0.234;
N = 30001;
x = x0 + linspace(-1.7343,1.99,N)*fwhm*2;
[y1,y2] = lorentzian(x,x0,fwhm,0);
[y3,y4] = lorentzian(x,x0,fwhm,0,pi/2);

% Compare
ok = areequal(y2,y3,1e-4,'abs');

