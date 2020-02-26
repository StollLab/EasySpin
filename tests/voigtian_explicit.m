function ok = test()

% Test voigtian against explitic construction

n = 2001;
x = linspace(-1,2,n);
x0 = 0.4;
fwhm = [0.1 0.07];

y1 = voigtian(x,x0,fwhm);
y2G = gaussian(x,x0,fwhm(1));
y2L = lorentzian(x,0.5,fwhm(2));

y2 = conv(y2G,y2L);

m = ceil(n/2);
y2 = y2(m:end-n+m);
y2 = y2*(x(2)-x(1));

ok = areequal(y1,y2,1e-4,'abs');
