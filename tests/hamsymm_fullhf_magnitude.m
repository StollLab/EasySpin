function ok = test()

% Test whether correct symmetry is detected for decreasing magnitude of
% hyperfine tensor.

% If sufficiently small, this test will fail. Failure point depends on
% threshold setting in hamsymm>eqeig. If tensor magnitude is too small,
% the numerically determined symmetry is higher than is actually is.

Sys.Nucs = '1H';
AFrame = [0 24 0]*pi/180;
Apv = [1 1 2];

R = erot(AFrame);
A = R*diag(Apv)*R.';

Ascale = [100 10 1 0.1 0.001];

for a = 1:numel(Ascale)

  Sys.A = A*Ascale(a);

  symm = hamsymm(Sys);
  ok(a) = strcmp(symm,'C2h');
end

