function ok = test()

% Rotation matrices with small beta (or beta near pi) obtained as products of
% several rotation matrices. These contain absolute rounding errors of the
% order of eps in R(3,1:2) and R(1:2,3), which must not be amplified by
% 1/sin(beta) in the rotation matrix reconstructed from the Euler angles.

rng(8812);

betaList = [1e-4 1e-6 1e-8 1e-10 1e-11];
nTrials = 20;

idx = 0;
for b = betaList
  for beta = [b pi-b]
    for t = 1:nTrials
      a1 = 2*pi*rand; b1 = pi*rand;
      a2 = 2*pi*rand; g1 = 2*pi*rand; g2 = 2*pi*rand;
      % Product equals erot(a2+g2,beta,g1)
      R = erot(0,0,g1)*erot(a1,b1,0)*erot(-a1,0,0)*erot(0,beta-b1,0)*erot(a2,0,g2);
      ang = eulang(R);
      idx = idx + 1;
      ok(idx) = norm(erot(ang)-R) < 1e-13 && abs(ang(2)-beta) < 1e-13;
    end
  end
end
