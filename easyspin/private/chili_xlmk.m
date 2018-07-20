% chili_xlmk   Computes the coefficients X^L_MK of the orienting
%              potential operator coefficients for the diffusion operator.
%
%     X = chili_xlk(Potential,R)
%
%          Potential.lambda   coefficients of potential expansion
%          Potential.L        L numbers of potential coefficients
%          Potential.M        M numbers of potential coefficients
%          Potential.K        K numbers of potential coefficients
%          R                  3 principal values of diffusion tensor

function X = chili_xlmk(Potential,R)

% Unpack potential parameters
idx = Potential.lambda~=0; % remove zero terms
lambda = Potential.lambda(idx);
Lpot = Potential.L(idx);
Mpot = Potential.M(idx);
Kpot = Potential.K(idx);

% Catch absence of potential
if isempty(lambda)
  X = cell(0,0);
  return
end

% Add terms needed to render potential real-valued
idx = Mpot~=0 | Kpot~=0;
lambda = [lambda conj(lambda(idx)).*(-1).^(Kpot(idx)-Mpot(idx))];
Lpot = [Lpot  Lpot(idx)];
Mpot = [Mpot -Mpot(idx)];
Kpot = [Kpot -Kpot(idx)];

maxLpot = max(Lpot);
maxLx = maxLpot*2;
nPotentialTerms = numel(lambda);

% Build lambda arrays
lam = cell(maxLpot+1,1);
for L = 0:maxLpot
  lam{L+1} = zeros(2*L+1);
end
for p = 1:nPotentialTerms
  lam{Lpot(p)+1}(Mpot(p)+Lpot(p)+1,Kpot(p)+Lpot(p)+1) = lambda(p);
end

% Get principal values of diffusion tensor
Rz = R(3);
Rp = (R(1)+R(2))/2;
Rd = (R(1)-R(2))/4;

% Initialize X output arrays
X = cell(maxLx+1,1);
for L = 0:maxLx
  X{L+1} = zeros(2*L+1);
end

% Define helper functions
Nm = @(L,K) sqrt((L+K-1)*(L+K)*(L-K+1)*(L-K+2));
Np = @(L,K) sqrt((L-K-1)*(L-K)*(L+K+1)*(L+K+2));
Mm = @(L,K) sqrt((L+K)*(L-K+1));
Mp = @(L,K) sqrt((L-K)*(L+K+1));

% Loop over LMK, calculate X^L_MK
for L = 0:maxLx
  for M = -L:L
    for K = -L:L
      
      % Calculate first term
      A = 0;
      if L<=maxLpot
        A = lam{L+1}(M+L+1,K+L+1) * (Rz*K^2 + Rp*(L*(L+1)-K^2));
        if Rd~=0
          if K+2<=L
            A = A + Rd*Nm(L,K+2)*lam{L+1}(M+L+1,K+2+L+1);
          end
          if K-2>=-L
            A = A + Rd*Np(L,K-2)*lam{L+1}(M+L+1,K-2+L+1);
          end
        end
      end
      
      % Calculate double sum (running only over non-zero lambdas)
      B = 0;
      for p1 = 1:nPotentialTerms
        L1 = Lpot(p1);
        M1 = Mpot(p1);
        K1 = Kpot(p1);
        lam1 = lambda(p1);
        if lam1==0, continue; end
        
        for p2 = 1:nPotentialTerms
          L2 = Lpot(p2);
          M2 = Mpot(p2);
          K2 = Kpot(p2);
          lam2 = lambda(p2);
          if lam2==0, continue; end
          
          % M-dependent 3j symbol: apply selection rules and calculate
          if M1-M+M2~=0, continue; end
          if abs(L1-L2)>L || L>L1+L2, continue; end
          if M1==0 && M2==0 && M==0
            if mod(L1+L2+L,2), continue; end
          end
          
          % Calculate all terms with K-dependent 3j symbols
          B_ = 0;
          if Rd~=0
            if K1-K+K2+2==0 && abs(K1+1)<=L1 && abs(K2+1)<=L2
              cpp = Mp(L1,K1)*Mp(L2,K2);
              if cpp~=0
                B_ = B_ + Rd*cpp*wigner3j(L1,L,L2,K1+1,-K,K2+1);
              end
            end
            if K1-K+K2-2==0 && abs(K1-1)<=L1 && abs(K2-1)<=L2
              cmm = Mm(L1,K1)*Mm(L2,K2);
              if cmm~=0
                B_ = B_ + Rd*cmm*wigner3j(L1,L,L2,K1-1,-K,K2-1);
              end
            end
          end
          if K1-K+K2==0
            if Rp~=0
              if abs(K1+1)<=L1 && abs(K2-1)<=L2
                cpm = Mp(L1,K1)*Mm(L2,K2);
                if cpm~=0
                  B_ = B_ + Rp*cpm*wigner3j(L1,L,L2,K1+1,-K,K2-1);
                end
              end
            end
            if Rz~=0 && K1~=0 && K2~=0
              B_ = B_ + Rz*K1*K2*wigner3j(L1,L,L2,K1,-K,K2);
            end
          end
          if B_==0, continue; end
                    
          B = B + lam1*lam2*wigner3j(L1,L,L2,M1,-M,M2)*B_;
          
        end
      end
      
      X{L+1}(M+L+1,K+L+1) = -1/2*A  - (2*L+1)/4*(-1)^(K-M)*B;
      
    end
  end
end
