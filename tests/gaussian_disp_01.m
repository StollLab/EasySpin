function [err,data] = test(opt,olddata)

% Check Gaussian dispersion explicit 1st deriv vs numerical deriv.
%=================================================================

% Line parameters
x0 = 2;
gamma = 1; % distance between inflection points
fwhm = gamma*sqrt(2*log(2));

x = linspace(-3,6,1001);

% Call gaussian() to get lineshapes
[dum,ydisp0] = gaussian(x,x0,fwhm,0);
[dum,ydisp1] = gaussian(x,x0,fwhm,1);

ydisp1_ = deriv(x,ydisp0);

if opt.Display
  shift = 0.0;
  plot(x,ydisp1,x,ydisp1_);
  legend('explicit','numerical deriv');
  legend boxoff
end

% Compare
err = ~areequal(ydisp1,ydisp1_,1e-3);

data =[];

