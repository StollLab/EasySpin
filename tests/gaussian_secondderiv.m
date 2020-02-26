function ok = test(opt)

% Check against explicit expression
%======================================================
x0 = 3; fwhm = 2;
x = linspace(-3,6,1001);

y = gaussian(x,x0,fwhm,2);

%gam = fwhm/sqrt(2*log(2));
%k = (x-x0)/gam;
%y0 = sqrt(2/pi)*4/gam^3*(4*k.^2-1).*exp(-2*k.^2);

gam = fwhm/sqrt(2*log(2))/2;
k = (x-x0)/gam;
y0 = sqrt(2/pi)*1/2/gam^3*(k.^2-1).*exp(-k.^2/2);

if opt.Display
  plot(x,y0,x,y,'r');
  legend('explicit','gaussian()');
  legend boxoff
end

ok = areequal(y,y0,1e-10,'rel');
