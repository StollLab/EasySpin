function ok = test()

% Rotation matrices generated with gamma = 0 have R(2,3) exactly zero and
% must give gamma exactly 0 (not a small negative number or a number close
% to 2*pi). For gamma = pi, R(2,3) is not exactly zero since sin(pi) is not
% exactly zero, so gamma is only required to be close to pi.

rng(7716);

nTrials = 50;

idx = 0;
for gamma0 = [0 pi]
  for t = 1:nTrials
    alpha0 = 2*pi*rand;
    beta0 = pi*rand;
    ang = eulang(erot(alpha0,beta0,gamma0));
    if gamma0==0
      gammaOk = ang(3)==0;
    else
      gammaOk = abs(ang(3)-pi)<1e-15;
    end
    idx = idx + 1;
    ok(idx) = gammaOk && areequal(ang(1:2),[alpha0 beta0],1e-10,'abs');
  end
end
