function [err,data] = test(opt,olddata)

% Check whether the Gaussian line width is correct
%-----------------------------------------------------
x0 = 0;
w = 0.1;
x = linspace(-1,1,10000)*2;
[yabs,ydisp] = gaussian(x,x0,w,0);

y_h = hilberttrans(yabs);

height1 = max(imag(y_h)) - min(imag(y_h));
height2 = max(ydisp)-min(ydisp);

% Compare
err = abs(height1-height2)>1e-2;

data =[];
