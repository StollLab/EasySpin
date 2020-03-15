function ok = test()

% Make sure Tkk = Nkk*Jp^k, with the appropriate normalization constant Nkk

J = 3/2;
Jp = sop(J,'+');

for k = 0:12
  Tkk = isto(J,[k k]);
  Nkk = (-1)^k/2^(k/2);
  Tkk_ref = Nkk*Jp^k;
  ok(k+1) = areequal(Tkk,Tkk_ref,1e-10,'abs');
end
