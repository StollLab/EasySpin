function ok = test(opt)

% Test whether the Lorentzian absorption lineshape and the derivative of
% its integral are consistent.

x = linspace(0,10,1e4);
x0 = 4;
fwhm = 1;

% Calculate Lorentzian absorption lineshape directly
y = lorentzian(x,x0,fwhm,0,pi/2);

% Calculate Lorentzian absorption lineshape via derivative of integral
y2 = lorentzian(x,x0,fwhm,-1,pi/2);
y2 = deriv(x,y2);

% Remove edge artifacts due to deriv()
y2([1 end]) = y([1 end]);

if opt.Display
  plot(x,y,x,y2);
  legend('direct','deriv of int');
end

ok = areequal(y,y2,1e-5,'rel');

end
