function [err,data] = test(opt,olddata)

% Check against explicit expression
%======================================================
x0 = 3; fwhm = 2;
x = linspace(-3,6,1001);

y = gaussian(x,x0,fwhm,1);

gam = fwhm/sqrt(2*log(2));
k = (x-x0)/gam;
y0 = -sqrt(2/pi)*4/gam^2*k.*exp(-2*k.^2);

if opt.Display
  plot(x,y0,x,y,'r');
  legend('explicit','gaussian()');
  legend boxoff
end

err = any(abs((y-y0)./y)>1e-10);

data =[];
