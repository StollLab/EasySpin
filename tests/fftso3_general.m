function ok = test()

Lmax = 5;
rng(623);

% Generate a general bandlimited function from Wigner expansion
idx = 0;
N = (Lmax+1)*(2*Lmax+1)*(2*Lmax+3)/3;
LMK = zeros(N,3);
vals = zeros(N,1);
for L = 0:Lmax
  for M = -L:L
    for K = -L:L
      idx = idx + 1;
      LMK(idx,:) = [L M K];
      vals(idx) = rand + 1i*rand;
    end
  end
end
f = @(a,b,c)fcn(LMK,vals,a,b,c);

% Calculate Fourier transform
[LMK_,vals_] = fftso3(f,Lmax);

ok = areequal(vals,vals_,1e-10,'abs') && areequal(LMK,LMK_);

end

function y = fcn(LMK,v,a,b,c)
y = 0;
for q = 1:numel(v)
  y = y + v(q)*wignerd(LMK(q,:),a,b,c);
end
end
