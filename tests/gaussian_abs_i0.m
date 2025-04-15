function ok = test(opt)

% Check Gaussian against derivative of its integral

x0 = 3; fwhm = 2;
x = linspace(-3,6,1001);

y0 = gaussian(x,x0,fwhm,-1);
y1 = gaussian(x,x0,fwhm,0);

y1_ = deriv(x,y0);

if opt.Display
  plot(x,y1,x,y1_,'r')
end

ok = areequal(y1,y1_,1e-3,'rel');
