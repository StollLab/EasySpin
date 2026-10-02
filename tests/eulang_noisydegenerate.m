function ok = test()

% Rotation matrices with beta = 0 or pi plus small noise (below the
% orthogonalization threshold) must be reconstructed to within the noise
% level. Noise in R(3,1:2) and R(1:2,3) makes alpha and gamma individually
% ill-defined, but alpha+gamma (beta=0) and gamma-alpha (beta=pi) must
% still be correct.

rng(5523);

noiseList = [1e-12 1e-11 1e-10 1e-9 1e-8];
beta0List = [0 pi];
nTrials = 20;

idx = 0;
for beta0 = beta0List
  for noise = noiseList
    for t = 1:nTrials
      R = erot(2*pi*rand,beta0,2*pi*rand) + noise*randn(3);
      ang = eulang(R);
      idx = idx + 1;
      ok(idx) = norm(erot(ang)-R) < 20*noise;
    end
  end
end
