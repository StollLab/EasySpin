function ok = test()

% Test Gaussian dispersion integral against explicit value

% gaussian() dispersion integral is only implemented using MATLAB's
% Symbolic Math Toolbox
if ~license('test','symbolic_toolbox')
  ok = true;
  return
end

x = 5;
x0 = 4;
fwhm = 1;

y = gaussian(x,x0,fwhm,-1,pi/2);

% Reference value
y_ref = 0.43443338893017619450;  % calculated with Mathematica
%{
Mathematica code:
x = 5; x0 = 4; fwhm = 1; sig = 
 fwhm/Sqrt[2*Log[2]]/2; k = (x - x0)/sig/Sqrt[2];
a = Sqrt[2/\[Pi]]*2/Sqrt[\[Pi]]/sig/2*sig*Sqrt[2]*k^2/2*
   HypergeometricPFQ[{1, 1}, {3/2, 2}, -k^2];
N[a, 20]
%}

ok = areequal(y,y_ref,1e-10,'rel');
