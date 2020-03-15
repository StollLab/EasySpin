function ok = test()

% Make sure operators returned from isto() satisfy Racah's commutation rule
% [Jm,Tkq] = sqrt((k+q)(k-q+1))*Tkq-1

Skq(1,:) = [5/2 3 2];
Skq(2,:) = [5 5 -1];

for i = 1:size(Skq,1)
  S = Skq(i,1);
  k = Skq(i,2);
  q = Skq(i,3);
  A = commute(sop(S,'-'),isto(S,[k q]));
  B = sqrt((k+q)*(k-q+1))*isto(S,[k q-1]);
  ok(i) = areequal(A,B,1e-12,'rel');
end
