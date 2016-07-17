function [err,data] = test(opt,olddata)

% Check Gaussian dispersion explicit 2nd deriv vs numerical deriv.
%=================================================================

% Line parameters
x0 = 2;
gamma = 1; % distance between inflection points
fwhm = gamma*sqrt(2*log(2));

x = linspace(-3,6,1001);

% Call gaussian() to get lineshapes
[dum,ydisp0] = gaussian(x,x0,fwhm,0);
[dum,ydisp2] = gaussian(x,x0,fwhm,2);

ydisp2_ = deriv(x,ydisp0);
ydisp2_ = deriv(x,ydisp2_);

if opt.Display
  shift = 0.0;
  plot(x,ydisp2,x,ydisp2_);
  legend('explicit','numerical deriv');
  legend boxoff
end

% Compare
err = ~areequal(ydisp2,ydisp2_,1e-2);

data =[];

