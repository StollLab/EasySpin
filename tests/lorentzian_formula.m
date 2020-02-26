function ok = test()

% Test: Check against formula (no derivative)

x0 = 3; fwhm = 2;
x = -3:6;
y = lorentzian(x,x0,fwhm,0);
gam = fwhm/sqrt(3);
xi = (x-x0)/gam;
y0 = 2/pi/sqrt(3)/gam./(1+4/3*xi.^2);

ok = areequal(y,y0,1e-10,'rel');

