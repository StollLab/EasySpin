function ok = test()

% Test with Wigner functions for all L, M, and K up to a maximum L.

A = 0.678+1.87i;

Lmax = 3;

idx = 0;
for L = 0:Lmax
  for M = -L:L
    for K = -L:L
      
      % Define single-D function, do FFT
      f = @(a,b,c) A*wignerd([L M K],a,b,c);
      c = fftso3(f,Lmax);
      
      % Extract expected non-zero element, and set to zeor
      A_ = c{L+1}(M+L+1,K+L+1);
      c{L+1}(M+L+1,K+L+1) = 0;
      
      idx = idx + 1;
      ok(idx) = areequal(A_,A,1e-10,'abs') && all(cellfun(@(x)all(x(:)==0),c));
    end
  end
end
