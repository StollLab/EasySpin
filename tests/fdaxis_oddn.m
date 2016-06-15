function [err,data] = test(opt,olddata)

% Test 2: Result test for odd N
%------------------------------------------------
dT = rand;

N = 345;
fx1 = fdaxis(dT,N);
fx2 = (-floor(N/2):floor(N/2))*(1/N/dT);

err = ~areequal(fx1,fx2);
data = [];
