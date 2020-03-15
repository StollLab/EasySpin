function ok = test()

% Make sure isto() operators satisfy [Jp,Tkq] = q*Tkq

S = 3/2;
k = 3;
Skq(:,3) = -k:k;
Skq(:,1) = S;
Skq(:,2) = k;

Jz = sop(S,'z');

for i = 1:size(Skq,1)
  S = Skq(i,1);
  k = Skq(i,2);
  q = Skq(i,3);
  Tkq = isto(S,[k q]);  
  A = commute(Jz,Tkq);
  B = q*Tkq;
  ok(i) = areequal(A,B,1e-10,'abs');
end
