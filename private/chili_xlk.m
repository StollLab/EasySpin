% chili_xlk    Computes the coefficients of the orienting
%              potential terms in the diffusion operator.
%
%     X = chili_xlk(Potential,R)
%
%          Potential.lambda   coefficients of potential expansion
%          Potential.L        L numbers of potential coefficients
%          Potential.K        K numbers of potential coefficients
%          R                  3 principal values of diffusion tensor

function X = chili_xlk(Potential,R)

% Implements Eq.(B6) from Lee et al. 1994
%   (here, we are using lambda instead of epsilon)
% Implements Eq.(A23) from Meirovitch et al. 1982
%   (missing prefactor 1/4 in front of sum)
% Implements Eq.(A7) from Earle et al. 1993
%   (missing - in front of first term in Eq. (A7))
%
% Missing minus in Earle: Is he correct and the others wrong?
% No, sign depends on definition of U = -k T ... or U = + k T...


% Shortcut if all lambda are zero
%---------------------------------------------------------------
if all(Potential.lambda==0) || isempty(Potential.lambda)
  logmsg(1,'Ordering potential: absent');
  X = zeros(0,0);
  return
end
logmsg(1,'Ordering potential: computing X(l,k) coefficients.');


% Get principal values of diffusion tensor
%---------------------------------------------------------------
Rz = R(3);
Rp = (R(1)+R(2))/2;
Rd = (R(1)-R(2))/4;


% Process lambda list
%---------------------------------------------------------------

maxL = max(Potential.L)*2;

if any(~isreal(Potential.lambda))
  error('Only real potential coefficients are allowed.');
end

%if any((KK<0)|(KK>LL))
%  error('L and K values of potential coefficients do not satisfy 0<=K<=L.');
%end

%if any(mod(LL,2) | mod(KK,2))
%  error('L and K of potential coefficients must be multiples of 2.');
%end

% Precompute 3j values of (L1,L,L2;0,0,0)
%---------------------------------------------------------------------
persistent jjj;
if (maxL>length(jjj)-1)
  jjj = zeros(maxL+1,maxL+1,maxL+1);
  for L = 0:maxL
    for L1 = 0:maxL
      for L2 = 0:maxL
        jjj((L1)+1,(L)+1,(L2)+1) = wigner3j(L1,L,L2,0,0,0);
      end
    end
  end
end

% Set up lambda matrix
%---------------------------------------------------------------------
% Layout of lambda, one L per row:
%  [ 0,  -,  -,  -,  -], L = 0
%  [-1,  0, +1,  -,  -], L = 1
%  [-2, -1,  0, +1, +2], L = 2
%  etc.
lambda = zeros(maxL+1,2*maxL+1);
for q = 1:numel(Potential.lambda)
  idxL = Potential.L(q)+1;
  lambda(idxL,+Potential.K(q)+idxL) = Potential.lambda(q);
  lambda(idxL,-Potential.K(q)+idxL) = Potential.lambda(q);
end

% Compute X matrix
%---------------------------------------------------------------------
% X^L_K = -1/2 * A  -  (2*L+1)/4 * B
% A: term linear in lambda
% B: term quadratic in lambda

Nm = @(L,K) sqrt((L+K-1)*(L+K)*(L-K+1)*(L-K+2));
Np = @(L,K) sqrt((L-K-1)*(L-K)*(L+K+1)*(L+K+2));
Mm = @(L,K) sqrt((L+K)*(L-K+1));
Mp = @(L,K) sqrt((L-K)*(L+K+1));

X = zeros(size(lambda));

for L = 0:maxL
  for K = -L:L
    
    A = lambda((L)+1,(K)+L+1) * (Rp*(L*(L+1)-K^2) + Rz*K^2);
    if (Rd~=0)
      if (K-2>=-L), A = A + lambda((L)+1,(K-2)+L+1)*Rd*Np(L,K-2); end
      if (K+2<=+L), A = A + lambda((L)+1,(K+2)+L+1)*Rd*Nm(L,K+2); end
    end

    B = 0;
    for L1 = 0:maxL
      for L2 = 0:maxL
        c0 = jjj((L1)+1,(L)+1,(L2)+1);
        if (c0==0), continue; end
        for K1 = -L1:L1
          if lambda((L1)+1,(K1)+L1+1)==0, continue; end
          for K2 = -L2:L2
            if lambda((L2)+1,(K2)+L2+1)==0, continue; end
            
            BB = 0;
            if (K1-K+K2==0)
              BB = Rz*K1*K2*wigner3j(L1,L,L2,K1,-K,K2);
              cpm = Mp(L1,K1)*Mm(L2,K2);
              if (cpm~=0)
                BB = BB + Rp*cpm*wigner3j(L1,L,L2,K1+1,-K,K2-1);
              end
            end
            
            if (Rd~=0)
              if (K1-K+K2+2==0)
                cpp = Mp(L1,K1)*Mp(L2,K2);
                if (cpp~=0)
                  BB = BB + Rd*cpp*wigner3j(L1,L,L2,K1+1,-K,K2+1);
                end
              end
              if (K1-K+K2-2==0)
                cmm = Mm(L1,K1)*Mm(L2,K2);
                if (cmm~=0)
                  BB = BB + Rd*cmm*wigner3j(L1,L,L2,K1-1,-K,K2-1);
                end
              end
            end
            
            B = B + BB * c0 * lambda((L1)+1,(K1)+L1+1) * lambda((L2)+1,(K2)+L2+1);
            
          end
        end
      end
    end

    value = -1/2*A  - (2*L+1)/4*(-1)^K*B;
    X((L)+1,(K)+L+1) = value;
    if abs(value)>0
    %fprintf('%d,%d: %f\n',L,K,value);
    end    
  end
end

return
