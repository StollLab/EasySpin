function [err,data] = test(opt,olddata)

% Test whether linefunction is normalized
%======================================================
x = linspace(-100,100,1e3); x0 = 34; w = 20;
y = gaussian(x,x0,w);
Area = trapz(x,y);

err = (abs(Area-1)>1e-10) | any(y<=0);
data = [];
