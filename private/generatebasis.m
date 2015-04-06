function basisList = generatebasis(Basis)

%%-------------------------------------------------------------------------
% Generate array of basis functions and their indices in full Liouvillian
%   
%   L  increasing from 0 to Lmax
%   M  decresing from L to -L
%   K  decreasing from L to 0
%   jK decreasing from +1 to -1 for nonzero K
%
%     L  M  K  jK
%%-------------------------------------------------------------------------

Leven = 0:2:Basis.evenLmax;
Lodd  = 1:2:Basis.oddLmax;
Llist = sort([Leven Lodd]);

Kmax = Basis.Kmax;
Mmax = Basis.Mmax;

iBasis = 1;
for L = Llist
  Mmx = min(L,Mmax);
  for M = Mmx:-1:-Mmx
    for K = min(L,Kmax):-1:0
      if (K~=0)
        basisList(iBasis,:)   = [L M K  1];
        basisList(iBasis+1,:) = [L M K -1];
        iBasis = iBasis + 2;
      else
        basisList(iBasis,:) = [L M K (-1)^L];
        iBasis = iBasis + 1;
      end
    end
  end
end

%{
% basis set ordering from S=1/2 code
jKmin = -1;
deltaK = 1;
iBasis = 1;
for L = Llist
  if (mod(L,2)==0), Lparity = +1; else Lparity = -1; end
  for jK = jKmin:2:1
    for K = 0:deltaK:min(L,Kmax)
      if ((K==0) && (Lparity~=jK)), continue; end
      Mmx = min(L,Mmax);
      for M = -Mmx:1:Mmx
        basisList(iBasis,:) = [L M K jK];
        iBasis = iBasis + 1;
      end % M
    end % K
  end % jK
end % L
%}

return
