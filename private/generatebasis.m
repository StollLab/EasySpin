function basisList = generatebasis(Basis,nSpin)

%%-------------------------------------------------------------------------
% Generate array of basis functions and their indices in full Liouvillian
%   
%   L  increasing from 0 to Lmax
%   M  decresing from L to -L
%   K  decreasing from L to 0
%   jK decreasing from +1 to -1 for nonzero K
%
%     L  M  K  jK  index
%%-------------------------------------------------------------------------

Leven = 0:2:Basis.evenLmax;
Lodd  = 1:2:Basis.oddLmax;
Llist = sort([Leven Lodd]);

Kmax = Basis.Kmax;
Mmax = Basis.Mmax;

iBasis = 1;
for L = Llist
  idxL = sum((2*(0:L-1)+1).^2)*nSpin;
  for M = min(L,Mmax):-1:max(-L,-Mmax)
    idxM = (L-M)*(2*L+1)*nSpin;
    for K = min(L,Kmax):-1:0
      idxK = 2*(L-K)*nSpin;
      idx = idxL + idxM + idxK + 1;
      if (K~=0)
        basisList(iBasis,:)   = [L M K  1 idx];
        basisList(iBasis+1,:) = [L M K -1 idx+nSpin];
        iBasis = iBasis + 2;
      else
        basisList(iBasis,:) = [L M K (-1)^L idx];
        iBasis = iBasis + 1;
      end
    end
  end
end

return
