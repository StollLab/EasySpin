function [err,data] = test(opt,olddata)

% vector of spin quantum numbers
qn = [1/2 3/2 1 1 1 1];
v = hsdim(qn);
err = ~areequal(v,prod(2*qn+1),0);
data = [];
