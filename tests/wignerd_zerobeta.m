function [err,data] = test(opt,olddata)

%===============================================================================
% Assert that wignerd gives identity matrix when alpha=beta=gamma=0
%===============================================================================

beta = 0;
ok = true;
for J = 0:0.5:6
  d = wignerd(J,beta);
  d0 = eye(2*J+1);
  ok = ok && areequal(d,d0,1e-16);
end
err = ~ok;

data = [];