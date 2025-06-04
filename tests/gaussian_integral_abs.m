function ok = test(opt)

% Test Gaussian absorption integral against explicit expression

x0 = 3;
fwhm = 2;
x = linspace(-3,6,1001);

y = gaussian(x,x0,fwhm,-1);

% Explicit expression
sig = fwhm/sqrt(2*log(2))/2;
k = (x-x0)/sig/sqrt(2);
y0 = (1+erf(k))/2;

if opt.Display
  plot(x,y0,x,y,'r');
  legend('explicit','gaussian()');
  legend boxoff
end

ok = areequal(y,y0,1e-10,'rel');
