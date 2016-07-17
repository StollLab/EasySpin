function [err,data] = test(opt,olddata)

% Check whether the Lorentzian line width is correct
%-----------------------------------------------------

fwhm = 1.2; % full width at half maximum

% Calculate Lorentzian lineshape
x0 = 0.234;
N = 30001;
x = x0 + linspace(-1.7343,1.99,N)*fwhm*1;
y = lorentzian(x,x0,fwhm,0);

% Determine fwhm from Lorentzian lineshape
y = y/max(y);
y = abs(y-0.5);
y1 = y; y1(x>x0) = 0.5;
[dum,idx1] = min(y1);
y2 = y; y2(x<x0) = 0.5;
[dum,idx2] = min(y2);
fwhm_est = x(idx2)-x(idx1);

% Compare
err = abs(fwhm_est-fwhm)>1e-4;

data =[];

