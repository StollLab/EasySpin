% jjjsymbol takes the maximum rank of the orientational basis as input and
% gives as output three 2-D arrays, storing the the following 3-j symbols:
%
%  rank 0: (L1 0 L2;  M1, -M1-M2, M2)
%  rank 1: (L1 1 L2;  M1, -M1-M2, M2)
%  rank 2: (L1 2 L2;  M1, -M1-M2, M2)
%
% L,M ordering (first to last): increasing L, and for each L, decreasing M.

function [jjj0,jjj1,jjj2] = jjjsymbol(evenLmax,oddLmax,computeRank1)

L = sort([0:2:evenLmax 1:2:oddLmax]);
nBasis = (max(L)+1)^2;

% Use cached results from previous runs if possible
%------------------------------------------------------------------
persistent cache
maxCache = 4;

if ~isempty(cache)
  for k = 1:numel(cache)
    c = cache(k);
    if c.evenLmax==evenLmax && ...
       c.oddLmax==oddLmax && ...
       c.computeRank1==computeRank1
      jjj0 = c.jjj0;
      jjj1 = c.jjj1;
      jjj2 = c.jjj2;
      return
    end
  end
end

% Rank-0 3j-symbols
%------------------------------------------------------------------
jjj0 = zeros(nBasis,nBasis);
for L1 = L
  for MK1 = L1:-1:-L1
    idx1 = L1^2 + (L1-MK1) + 1;
    for L2 = L
      for MK2 = L2:-1:-L2
        idx2 = L2^2 + (L2-MK2) + 1;
        
        if (L1 ~= L2) && (MK1 ~= -MK2), continue; end
        
        jjj0(idx1,idx2) = wigner3j(L1,0,L2,MK1,0,MK2);
        
      end
    end
  end
end

% Rank-1 3j-symbols
%------------------------------------------------------------------
if computeRank1
  jjj1 = zeros(nBasis,nBasis);
  for L1 = L
    for MK1 = L1:-1:-L1
      idx1 = L1^2 + (L1-MK1) + 1;
      for L2 = L
        for MK2 = L2:-1:-L2
          idx2 = L2^2 + (L2-MK2) + 1;
          
          if abs(-MK1-MK2) > 1, continue; end
          if abs(L1-L2) > 1, continue; end
          
          jjj1(idx1,idx2) = wigner3j(L1,1,L2,MK1,-MK2-MK1,MK2);
          
        end
      end
    end
  end
else
  jjj1 = [];
end

% Rank-2 3j-symbols
%------------------------------------------------------------------
jjj2 = zeros(nBasis,nBasis);
for L1 = L
  for MK1 = L1:-1:-L1
    idx1 = L1^2 + (L1-MK1) + 1;
    for L2 = L
      for MK2 = L2:-1:-L2
        idx2 = L2^2 + (L2-MK2) + 1;
        
        if abs(-MK1-MK2) > 2, continue; end
        if abs(L1-L2) > 2, continue; end
        
        jjj2(idx1,idx2) = wigner3j(L1,2,L2,MK1,-MK2-MK1,MK2);
        
      end
    end
  end
end

% Store in internal cache
%------------------------------------------------------------------
% If cache size is going to be exceeded, drop oldest cached result.
if numel(cache)>=maxCache
  cache(1) = [];
end

% Store result in internal cache.
idx = numel(cache)+1;
cache(idx).evenLmax = evenLmax;
cache(idx).oddLmax = oddLmax;
cache(idx).computeRank1 = computeRank1;
cache(idx).jjj0 = jjj0;
cache(idx).jjj1 = jjj1;
cache(idx).jjj2 = jjj2;

return
