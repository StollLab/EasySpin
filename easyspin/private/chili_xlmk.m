% chili_xlmk   Compute the coefficients X^L_MK of the potential-dependent
%              part of the diffusion operator.
%
%     X = chili_xlmk(Potential,R)
%
% Input:
%   Potential.lambda   coefficients of potential expansion
%   Potential.L        L numbers of potential coefficients
%   Potential.M        M numbers of potential coefficients
%   Potential.K        K numbers of potential coefficients
%   R                  3 principal values of diffusion tensor
%
% Output:
%   X                  cell array of X^L_MK

function X = chili_xlmk(Potential,R)

% Check for absence of potential
%-------------------------------------------------------------------------------
if isempty(Potential.lambda) || all(Potential.lambda==0)
  X = cell(0,0);
  return
end

% Unpack potential coefficients, remove zero coefficients
%-------------------------------------------------------------------------------
Lpot = Potential.L(:);
Mpot = Potential.M(:);
Kpot = Potential.K(:);
lambda = Potential.lambda(:);
rmv = lambda==0;
if any(rmv)
  Lpot(rmv) = [];
  Mpot(rmv) = [];
  Kpot(rmv) = [];
  lambda(rmv) = [];
end

% Add terms needed to render potential real-valued
%-------------------------------------------------------------------------------
idx = Mpot~=0 | Kpot~=0;
lambda = [lambda; conj(lambda(idx)).*(-1).^(Kpot(idx)-Mpot(idx))];
Lpot = [Lpot;  Lpot(idx)];
Mpot = [Mpot; -Mpot(idx)];
Kpot = [Kpot; -Kpot(idx)];

nPotentialTerms = numel(lambda);
maxLpot = max(Lpot);
allMzero = all(Mpot==0);

% Set up lambda arrays
%-------------------------------------------------------------------------------
lam_ = cell(maxLpot+1,1);
for L = 0:maxLpot
  lam_{L+1} = zeros(2*L+1);
end
for p = 1:nPotentialTerms
  L1_ = Lpot(p)+1;
  lam_{L1_}(Mpot(p)+L1_,Kpot(p)+L1_) = lambda(p);
end
lam = @(L,M,K) lam_{L+1}(M+L+1,K+L+1);

% Get principal values of diffusion tensor
%-------------------------------------------------------------------------------
switch numel(R)
  case 1 % isotropic
    R = [R R R];
  case 3 % three principal values
  otherwise
    error('Provide 3 principal values of diffusion tensor as second input.');
end
Rz = R(3);
Rp = (R(1)+R(2))/2;
Rd = (R(1)-R(2))/4;

% Initialize X output array
%-------------------------------------------------------------------------------
maxLx = maxLpot*2;
X = cell(maxLx+1,1);
for L = 0:maxLx
  X{L+1} = zeros(2*L+1);
end

% Calculate X^L_MK coefficients
%-------------------------------------------------------------------------------
% X^L_MK = -1/2 * A  -  1/4 * (2*L+1) * (-1)^(K-M) * B
%   A: three terms linear in lambda
%   B: double sum of terms quadratic in lambda

% Define helper functions
cm = @(L,K) sqrt(L*(L+1)-K*(K-1));
cp = @(L,K) sqrt(L*(L+1)-K*(K+1));

% Loop over LMK, calculate X^L_MK
for L = 0:maxLx
  if allMzero, Mrange = 0; else, Mrange = -L:L; end
  for M = Mrange
    for K = -L:L
      
      % Calculate A terms
      A = 0;
      if L<=maxLpot
        if Rd~=0
          if K+2<=L
            A = A + Rd*cm(L,K+1)*cm(L,K+2)*lam(L,M,K+2);
          end
          if K-2>=-L
            A = A + Rd*cp(L,K-1)*cp(L,K-2)*lam(L,M,K-2);
          end
        end
        A = A + lam(L,M,K) * (Rp*(L*(L+1)-K^2) + Rz*K^2);
      end
      
      % Calculate B terms
      % (double sum; running only over non-zero potential coefficients)
      B = 0;
      for p1 = 1:nPotentialTerms
        lam1 = lambda(p1);
        if lam1==0, continue; end
        L1 = Lpot(p1);
        M1 = Mpot(p1);
        K1 = Kpot(p1);
        
        for p2 = 1:nPotentialTerms
          lam2 = lambda(p2);
          if lam2==0, continue; end
          L2 = Lpot(p2);
          M2 = Mpot(p2);
          K2 = Kpot(p2);
          
          % M-dependent 3j symbol: apply selection rules
          if M1-M+M2~=0, continue; end
          if abs(L1-L2)>L || L>L1+L2, continue; end
          if M1==0 && M2==0 && M==0
            if mod(L1+L2+L,2), continue; end
          end
          
          % Calculate all terms with K-dependent 3j symbols
          B_ = 0;
          if Rd~=0
            if K1-K+K2+2==0 && abs(K1+1)<=L1 && abs(K2+1)<=L2
              cpp = cp(L1,K1)*cp(L2,K2);
              if cpp~=0
                B_ = B_ + Rd*cpp*wigner3j(L1,L,L2,K1+1,-K,K2+1);
              end
            end
            if K1-K+K2-2==0 && abs(K1-1)<=L1 && abs(K2-1)<=L2
              cmm = cm(L1,K1)*cm(L2,K2);
              if cmm~=0
                B_ = B_ + Rd*cmm*wigner3j(L1,L,L2,K1-1,-K,K2-1);
              end
            end
          end
          if K1-K+K2==0
            if Rz~=0 && K1~=0 && K2~=0
              B_ = B_ + Rz*K1*K2*wigner3j(L1,L,L2,K1,-K,K2);
            end
            if Rp~=0
              if abs(K1+1)<=L1 && abs(K2-1)<=L2
                cpm = cp(L1,K1)*cm(L2,K2);
                if cpm~=0
                  B_ = B_ + Rp*cpm*wigner3j(L1,L,L2,K1+1,-K,K2-1);
                end
              end
            end
          end
          if B_==0, continue; end
          
          B = B + lam1*lam2*wigner3j(L1,L,L2,M1,-M,M2)*B_;
          
        end
      end
      
      % Combine A and B terms, and store
      X_ = -1/2*A  - (2*L+1)/4*(-1)^(K-M)*B;
      X{L+1}(M+L+1,K+L+1) = X_;
      
    end
  end
end
