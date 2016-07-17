function [err,data] = test(opt,olddata)

% Test 3: explicit formula for O_6^1 and many S values
%============================================================
kk = [6 4]; qq = [-1 -3];
for ii = 1:numel(kk)
  k = kk(ii);
  q = qq(ii);
  for S = 3:0.5:5
    s = S*(S+1);
    [Sp,Sm,Sz,Se] = sop(S,'+','-','z','e');
    if (k==6)&&(q==-1)
      a = 33*Sz^5 - (30*s-15)*Sz^3 + (5*s^2-10*s+12)*Sz;
    elseif (k==4)&&(q==-3)
      a = Sz;
    end
    b = Sp^abs(q) - Sm^abs(q);
    if q>0, c = 1/2; else c = 1/2i; end
    Op0 = c*(a*b + b*a)/2;
    Op1 = stev(S,k,q);
    err = any(abs(Op0(:)-Op1(:))>1e-6);
    if (err), break; end
  end
end
data = [];
