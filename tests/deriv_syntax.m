function [err,data] = test(opt,olddata)

% Test 1: syntax
%================================================================
N = 654;
y = rand(1,N);
x = 1:N;
a = deriv(y);
b = deriv(x,y);
err = ~areequal(a,b);
data = [];
