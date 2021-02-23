function ok = test()

% Difficult test case: R numerically very close to identity
% This should give all-zero angles.
% Specifically, alpha and gamma should be close to zero, and not close to 2*pi.

gv = [
   0.00006080549145  -0.08505970942061  -0.99637585385033
  -0.00049522568221   0.99637573095109  -0.08505972915078
   0.99999987552710   0.00049860301059   0.00001846136337];

R = gv*gv';
angles = eulang(R);

ok = areequal(angles,[0 0 0],1e-12,'abs');
