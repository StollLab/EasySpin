function [err,data] = test(opt,olddata)

%===============================================================================
% Assert that wignerd gives identity matrix when alpha=beta=gamma=0
%===============================================================================

J = 4;
beta = 0;

d = wignerd(J,beta);
d0 = eye(2*J+1);

err = ~areequal(d,d0,1e-16);

data = [];