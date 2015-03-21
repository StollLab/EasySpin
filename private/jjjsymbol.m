% jjjsymbol takes the maximum rank of the orientational basis as input and
% gives as output two 2-D arrays, one storing the 3-j symbols for RBO rank
% l=0 and one for l=2.

function [jjj0,jjj1,jjj2] = jjjsymbol(LLKM)

% construct array containing indexing all combinations of L and M (or K)
%   [L M]

% Basis ordering (first to last): increasing L, and for each
% L, decreasing M.

Lmax = max(LLKM(1),LLKM(2));
nBasis = (Lmax+1)^2;

% calculate the 3j-symbols
%
%  (L1 0 L2;  M1, -M1-M2, M2)
%  (L1 2 L2;  M1, -M1-M2, M2)


% rank-0 3j-symbols

jjj0 = zeros(nBasis,nBasis);
for L1 = 0:Lmax
  for MK1 = L1:-1:-L1
    idx1 = L1^2 + (L1-MK1) + 1;
    for L2 = 0:Lmax
      for MK2 = L2:-1:-L2
        idx2 = L2^2 + (L2-MK2) + 1;
        
        if (L1 ~= L2) && (MK1 ~= -MK2), continue; end;
        
        jjj0(idx1,idx2) = wigner3j(L1,0,L2,MK1,0,MK2);
        
      end
    end
  end
end


% rank-1 3j-symbols

jjj1 = zeros(nBasis,nBasis);
for L1 = 0:Lmax
  for MK1 = L1:-1:-L1
    idx1 = L1^2 + (L1-MK1) + 1;
    for L2 = 0:Lmax
      for MK2 = L2:-1:-L2
        idx2 = L2^2 + (L2-MK2) + 1;
        
        if abs(-MK1-MK2) > 1, continue; end;
        if abs(L1-L2) > 1, continue; end;
        
        jjj1(idx1,idx2) = wigner3j(L1,1,L2,MK1,-MK2-MK1,MK2);
        
      end
    end
  end
end
    

% rank-2 3j-symbols

jjj2 = zeros(nBasis,nBasis);
for L1 = 0:Lmax
  for MK1 = L1:-1:-L1
    idx1 = L1^2 + (L1-MK1) + 1;
    for L2 = 0:Lmax
      for MK2 = L2:-1:-L2
        idx2 = L2^2 + (L2-MK2) + 1;
        
        if abs(-MK1-MK2) > 2, continue; end;
        if abs(L1-L2) > 2, continue; end;
        
        jjj2(idx1,idx2) = wigner3j(L1,2,L2,MK1,-MK2-MK1,MK2);
        
      end
    end
  end
end


return
