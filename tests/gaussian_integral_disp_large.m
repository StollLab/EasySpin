function ok = test()

% Test Gaussian dispersion integral against explicit value (large k, asymptotic expansion)

x = 10;
x0 = 4;
fwhm = 1;

y = gaussian(x,x0,fwhm,-1,pi/2);

% Reference value
y_ref = 1.0443369788200551275;  % calculated with MATLAB Symbolic Math Toolbox
%{
MATLAB code:
x = sym(10); x0 = sym(4); fwhm = sym(1);
sig = fwhm/sqrt(2*log(sym(2)))/2; k = (x-x0)/sig/sqrt(sym(2));
a = sqrt(2/sym(pi))*2/sqrt(sym(pi))/sig/2*sig*sqrt(sym(2))*k^2/2*...
  hypergeom([1,1],[sym(3)/2,2],-k^2);
vpa(a,20)
%}

ok = areequal(y,y_ref,1e-10,'rel');
