function [err,data] = test(opt,olddata)

% Check whether the Gaussian line width is correct
%-----------------------------------------------------

fwhm = 1.1; % full width at half maximum

% Calculate Gaussian lineshape
x0 = 0.234;
N = 30001;
x = x0 + linspace(-1.7343,1.99,N)*fwhm*2;
[y1,y2] = gaussian(x,x0,fwhm,0);
[y3,y4] = gaussian(x,x0,fwhm,0,pi/2);

% Compare
err = abs(y2-y3)>1e-4;

data =[];

