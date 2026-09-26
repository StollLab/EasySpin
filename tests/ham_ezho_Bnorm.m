function ok = test()

% Assert that the normalization of the field-dependent tensor is correct for
% orders lB=1..8. For B along z and a single term with (lB,lS,l)=(lB,lB,0),
% H = b^lB*(-1)^lB/sqrt(alpha*(2*lB+1))*T_lB0(S), with alpha = (2*lB-1)!!/lB!.

S = 4;
b = 1.2;
Sp = sop(S,'+');
Sm = sop(S,'-');
threshold = 1e-8;

for lB = 8:-1:1
  % T_lB0(S), by lowering from T_lB,lB = (-1)^lB*2^(-lB/2)*S+^lB
  T = (-1)^lB*2^(-lB/2)*Sp^lB;
  for m = lB:-1:1
    T = (Sm*T-T*Sm)/sqrt((lB+m)*(lB-m+1));
  end
  alpha = prod(1:2:2*lB-1)/factorial(lB);
  Href = b^lB*(-1)^lB/sqrt(alpha*(2*lB+1))*T;
  
  Sys = struct('S',S);
  Sys.(sprintf('Ham%i%i0',lB,lB)) = 1;
  H = ham_ezho(Sys,[0 0 b]);
  
  ok(lB) = areequal(H,full(Href),threshold,'abs');
end
