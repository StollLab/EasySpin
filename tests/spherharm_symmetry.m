function ok = test(opt,olddata)

% Test symmetry (L,M) <-> (L,-M)
%-------------------------------------------------------------------------------

theta = 0.2134;
phi  = 0.8765;

idx = 0;
for L = 0:5
  for M = 0:L
    Yp = spherharm(L,+M,theta,phi);
    Ym = spherharm(L,-M,theta,phi);
    idx = idx + 1;
    ok(idx) = areequal(Ym,(-1)^M*conj(Yp),1e-10,'rel');
  end
end
