function [err,data] = test(opt,olddata)

% Test 3: Result test for even N
%------------------------------------------------
dT = 10;

N = 532;
fx1 = fdaxis(dT,N);
fx2 = linspace(-N/2,N/2-1,N)*(1/N/dT);

err = ~areequal(fx1,fx2);
data = [];
