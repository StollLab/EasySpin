function ok = test()

% Small but non-zero beta, and beta just below pi, must not be treated as
% degenerate. The rotation matrix must be recovered exactly.

alpha = 0.7;
gamma = 0.3;
betalist = [1e-4 1e-6 1e-9 pi-1e-4 pi-1e-6 pi-1e-9];

for k = numel(betalist):-1:1
  R = erot(alpha,betalist(k),gamma);
  ang = eulang(R);
  ok(k) = areequal(erot(ang),R,1e-12,'abs') && ...
          areequal(ang,[alpha betalist(k) gamma],1e-6,'abs');
end
