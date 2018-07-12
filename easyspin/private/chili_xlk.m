% chili_xlk    Computes the coefficients X^L_0K of the orienting
%              potential operator coefficients for the diffusion operator.
%
%     X = chili_xlk(Potential,R)
%
%          Potential.lambda   coefficients of potential expansion
%          Potential.L        L numbers of potential coefficients
%          Potential.K        K numbers of potential coefficients
%          R                  3 principal values of diffusion tensor
%
%    chili_xlk does not support potentials with any Potential.M ~=0
%    In that case, use chili_xlmk.

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

if any(Potential.M)
  error('chili_xlk does not support potential coefficients with M~=0.');
end


% Shortcut if all lambda are zero
%---------------------------------------------------------------
if ~isfield(Potential,'lambda') || ~any(Potential.lambda)
  X = zeros(0,0);
  return
end


% Get principal values of diffusion tensor
%---------------------------------------------------------------
Rz = R(3);
Rp = (R(1)+R(2))/2;
Rd = (R(1)-R(2))/4;


% Process lambda list
%---------------------------------------------------------------

maxLpot = max(Potential.L(Potential.lambda~=0));
maxLx = maxLpot*2;

if any(~isreal(Potential.lambda))
  error('Only real potential coefficients are allowed.');
end

% Set up lambda matrix
%---------------------------------------------------------------------
% Layout of lambda, one L per row:
%  [ 0,  -,  -,  -,  -], L = 0
%  [-1,  0, +1,  -,  -], L = 1
%  [-2, -1,  0, +1, +2], L = 2
%  etc.
lambda = zeros(maxLpot+1,2*maxLpot+1);
for q = 1:numel(Potential.lambda)
  idxL = Potential.L(q)+1;
  lambda(idxL,+Potential.K(q)+idxL) = Potential.lambda(q);
  % assure that potential is real-valued
  lambda(idxL,-Potential.K(q)+idxL) = Potential.lambda(q)*(-1)^Potential.K(q);
end

% Compute X^L_K factors
%---------------------------------------------------------------------
% X^L_K = -1/2 * A  -  (2*L+1)/4 * B
% A: term linear in lambda
% B: term quadratic in lambda

% Define helper functions
Nm = @(L,K) sqrt((L+K-1)*(L+K)*(L-K+1)*(L-K+2));
Np = @(L,K) sqrt((L-K-1)*(L-K)*(L+K+1)*(L+K+2));
Mm = @(L,K) sqrt((L+K)*(L-K+1));
Mp = @(L,K) sqrt((L-K)*(L+K+1));

X = zeros(maxLx+1,2*maxLx+1);

for L = 0:maxLx
  for K = -L:L
    
    % Calculate first term
    A = 0;
    if L<=maxLpot
      A = lambda(L+1,K+L+1) * (Rp*(L*(L+1)-K^2) + Rz*K^2);
      if Rd~=0
        if K+2<=L
          A = A + Rd*Nm(L,K+2)*lambda(L+1,K+2+L+1);
        end
        if K-2>=-L
          A = A + Rd*Np(L,K-2)*lambda(L+1,K-2+L+1);
        end
      end
    end
    
    % Calculate second term
    B = 0;
    for L1 = 0:maxLpot
      for L2 = 0:maxLpot
        
        if abs(L1-L2)>L || L>L1+L2, continue; end % 3j selection rule
        if mod(L1+L2+L,2), continue; end % 3j selection rule
        jjj0 = wigner3j(L1,L,L2,0,0,0);
        if jjj0==0, continue; end
        
        for K1 = -L1:L1
          lam1 = lambda(L1+1,K1+L1+1);
          if lam1==0, continue; end
          for K2 = -L2:L2
            lam2 = lambda(L2+1,K2+L2+1);
            if lam2==0, continue; end
            
            B_ = 0;
            if K1-K+K2==0 % 3j selection rule
              B_ = Rz*K1*K2*wigner3j(L1,L,L2,K1,-K,K2);
              cpm = Mp(L1,K1)*Mm(L2,K2);
              if (cpm~=0)
                B_ = B_ + Rp*cpm*wigner3j(L1,L,L2,K1+1,-K,K2-1);
              end
            end
            
            if Rd~=0
              if K1-K+K2+2==0  % 3j selection rule
                cpp = Mp(L1,K1)*Mp(L2,K2);
                if cpp~=0
                  B_ = B_ + Rd*cpp*wigner3j(L1,L,L2,K1+1,-K,K2+1);
                end
              end
              if K1-K+K2-2==0  % 3j selection rule
                cmm = Mm(L1,K1)*Mm(L2,K2);
                if cmm~=0
                  B_ = B_ + Rd*cmm*wigner3j(L1,L,L2,K1-1,-K,K2-1);
                end
              end
            end
            if B_==0, continue; end
            
            B = B + lam1*lam2*jjj0*B_;
            
          end
        end
      end
    end
    
    X(L+1,K+L+1) = -1/2*A  - (2*L+1)/4*(-1)^K*B;
  end
end

return
