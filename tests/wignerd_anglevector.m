function [err,data] = test(opt,olddata)

%============================================================================
% Check whether single-matrix-element calls to wignerd() work with arrays
% of angles.
%============================================================================

rng_(62323,'twister');
n1 = 3;
n2 = 10;

alpha = rand(n1,n2)*2*pi;
beta = rand(n1,n2)*pi;
gamma = rand(n1,n2)*2*pi;

ok = [];
idx = 0;
for J = 0:0.5:5
  for m1 = -J:J
    for m2 = -J:J
      D = wignerd([J m1 m2],alpha,beta,gamma);
      idx = idx + 1;
      ok(idx) = all(size(D)==[n1 n2]);
    end
  end
end
err = any(~ok);

data = [];