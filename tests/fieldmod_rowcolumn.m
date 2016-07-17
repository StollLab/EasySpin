function [err,data] = test(opt,olddata)

%=======================================================
% Make sure fieldmod takes row and column vectors
%=======================================================
x = linspace(-10,10,1e3);
y = gaussian(x,-2,1);
dx = x(2)-x(1);
y1 = fieldmod(x,y,dx);
y2 = fieldmod(x,y.',dx);
y3 = fieldmod(x.',y,dx);
y4 = fieldmod(x.',y.',dx);

err = 0;
data = [];
