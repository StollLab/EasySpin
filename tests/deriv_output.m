function [err,data] = test(opt,olddata)

% Test 2: output
%================================================================
N = 2345;
y = rand(1,N);
x = 1:N;
a = deriv(x,y);
err = numel(a)~=numel(y);
data = [];
