function ok = test()

% TExplicit formula for O_6^1 and many S values

kq = [6 1; 4 3];
for ii = 1:size(kq,1)
  k = kq(ii,1);
  q = kq(ii,2);
  for S = 3:0.5:7
    s = S*(S+1);
    [Sp,Sm,Sz,Se] = sop(S,'+','-','z','e');
    if k==6 && q==1
      a = 33*Sz^5 - (30*s-15)*Sz^3 + (5*s^2-10*s+12)*Sz;
    elseif k==4 && q==3
      a = Sz;
    end
    b = Sp^q + Sm^q;
    Op0 = (a*b + b*a)/4;
    Op1 = stev(S,[k,q]);
    ok = areequal(Op0,Op1,1e-6,'abs');
    if ~ok, break; end
  end
end
