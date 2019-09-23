function [err,data] = test(opt,olddata)

% Test 1: Interface test
%------------------------------------------------
dT = 0.45;
N = 200;
tx = (0:N-1)*dT;
fx1 = fdaxis(dT,N);
fx2 = fdaxis(tx);
fx3 = fdaxis(tx + rand);
err = ~areequal(fx1,fx2,1e-10,'abs');
data = [];
