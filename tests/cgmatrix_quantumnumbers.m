function ok = test()

% Assert that cgmatrix() returns correct lists of quantum numbers for
% the coupled and the uncoupled basis.

S1 = 1;
S2 = 3/2;
[~,Sm,m12] = cgmatrix(S1,S2);

nStates = (2*S1+1)*(2*S2+1);

% Construct reference list of coupled quantum numbers
Sm_ref = zeros(nStates,2);
ic = 1;
for Stot = S1+S2:-1:abs(S1-S2)
  for mStot = Stot:-1:-Stot
    Sm_ref(ic,:) = [Stot mStot];
    ic = ic + 1;
  end
end

% Construct reference list of uncoupled quantum numbers
m12_ref = zeros(nStates,2);
iu = 1;
for m1 = S1:-1:-S1
  for m2 = S2:-1:-S2
    m12_ref(iu,:) = [m1 m2];
    iu = iu + 1;
  end
end

ok(1) = areequal(Sm,Sm_ref,1e-16,'abs');
ok(2) = areequal(m12,m12_ref,1e-16,'abs');
