function [err,data] = test(opt,olddata)

% Test 4: correct values
%================================================================
N = 1000;
x = 1:N; dx = x(2) - x(1);
y = rand(1,N);
dydx = diff(y)/dx;
dydx = (dydx([1 1:end]) + dydx([1:end end]))/2;
a = deriv(x,y);
err = ~areequal(dydx,a);
data = [];
