function ok = test()

% Test with Wigner functions for all L, M, and K up to a maximum L.

A = 0.678+1.87i; % arbitrary amplitude

Lmax = 4;

idx = 0;
for L = 0:Lmax
  for M = -L:L
    for K = -L:L
      
      f = @(a,b,c) A*wignerd([L M K],a,b,c);
      
      [LMK,A_] = fftso3(f,Lmax);
      
      idx = idx + 1;
      ok(idx) = all(LMK==[L M K]) && areequal(A_,A,1e-10,'abs');
    end
  end
end
