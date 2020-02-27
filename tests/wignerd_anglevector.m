function ok = test()

% Check whether single-matrix-element calls to wignerd() work with arrays
% of angles.

rng(62323);
n1 = 3;
n2 = 10;

alpha = rand(n1,n2)*2*pi;
beta = rand(n1,n2)*pi;
gamma = rand(n1,n2)*2*pi;

idx = 1;
for J = 0:0.5:5
  ok(idx) = true;
  for m1 = -J:J
    for m2 = -J:J
      D = wignerd([J m1 m2],alpha,beta,gamma);
      ok(idx) = ok(idx) && all(size(D)==[n1 n2]);
    end
  end
  idx = idx + 1;
end
