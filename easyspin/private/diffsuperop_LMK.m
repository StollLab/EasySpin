% This function computes the diffusion superoperator matrix in the LMK basis,
% for a diagonal diffusion tensor. It is used by chili().
%
% It does not require any particular order of orientational basis functions.
%
% Syntax:
%   Gamma = diffsuperop_LMK(basis,R,XLK,Method)
%
% Input:
%   basis     structure containing the lists for the quantum numbers for
%             the orientational basis in basis.L, basis.M, and basis.K
%   R         array with 3 principal values of diffusion tensor (s^-1),
%             or single value if diffusion is isotropic
%   Potential (optional) structure containing information on potential
%             .X        orientational potential coefficients, as returned by
%                       chili_xlmk (assumed all zero if not given)
%             .lambda   potential coefficients
%             .L,.M.,K  corresponding function indices
%   XLMK      diffusion operator coefficients, as provided by chili_xlmk
%   Method    (optional)
%             1: method for diagonal diffusion tensors
%             2: method for diagonal diffusion tensors (faster than 1) (default)
%
% Output:
%   Gamma    diffusion superoperator matrix in the LMK basis (s^-1), sparse

function Gamma = diffsuperop_LMK(basis,R,XLMK,Potential,Method)

if isfield(basis,'jK') && ~isempty(basis.jK)
  error('This function expects an LMK basis, without jK.');
end

if nargin<4, Method = []; end

usePotential = nargin>2 && any(Potential.lambda~=0);
if usePotential
  % Create a list of non-zero X coefficients and their indices Lx, Mx, Kx
  % from Potential.X
  Lx = [];
  Mx = [];
  Kx = [];
  X = [];
  for L_ = 0:numel(XLMK)-1
    [r,c,val] = find(XLMK{L_+1});
    Mx = [Mx r-L_-1];
    Kx = [Kx c-L_-1];
    Lx = [Lx L_*ones(1,numel(c))];
    X = [X val];
  end
end

if isempty(Method), Method = 1; end

switch numel(R)
  case 1 % isotropic value
    Rd = 0;
    Rperp = R;
    Rz = R;
    isAxial = Rd==0;
  case 3 % 3 principal values
    if isempty(Method), Method = 2; end
    Rd = (R(1)-R(2))/4;
    Rperp = (R(1)+R(2))/2;
    Rz = R(3);
    isAxial = Rd==0;
  otherwise
    error('Wrong number of elements in diffusion tensor.');
end

% Calculate diffusion operator
%-------------------------------------------------------------------------------
L = basis.L;
M = basis.M;
K = basis.K;
nBasis = numel(L);

if Method==1
  % Method 1: Calculate the diffusion operator directly (fast).
  %-----------------------------------------------------------------------------
  
  % Treat the cases of isotropic and axial diffusion tensors in the absence of
  % an orientational potential. In these cases, the diffusion operator matrix is
  % diagonal.
  if isAxial && ~usePotential
    diagonal = Rperp*(L.*(L+1)-K.^2) + Rz*K.^2;
    Gamma = spdiags(diagonal,0,nBasis,nBasis);
    return
  end
  
  % Only diagonal and upper triangular part of matrix are explicitly calulcated,
  % since Gamma is symmetric if R is diagonal.
  idx = 0;
  for b1 = 1:nBasis
    L1 = L(b1);
    M1 = M(b1);
    K1 = K(b1);
    
    for b2 = b1:nBasis
      M2 = M(b2);
      if M1~=M2, continue; end
      L2 = L(b2);
      K2 = K(b2);
      
      % Potential-independent part
      val_ = 0;
      if L1==L2
        LL2 = L2*(L2+1);
        if K1==K2
          val_ = Rperp * (LL2-K1^2) + Rz * K1^2;
        end
        if Rd~=0
          if K1==K2+2
            val_ = Rd * sqrt((LL2-K2*(K2+1))*(LL2-(K2+1)*(K2+2)));
          elseif K1==K2-2
            val_ = Rd * sqrt((LL2-K2*(K2-1))*(LL2-(K2-1)*(K2-2)));
          end
        end
      end
      
      % Potential-dependent part
      if usePotential
        valp_ = 0;
        for p = 1:numel(X)
          if Kx(p)~=K1-K2, continue; end
          if Mx(p)~=M1-M2, continue; end
          if abs(L1-Lx(p))>L2 || L2>L1+Lx(p), continue; end
          
          % Calculate M- and K-dependent 3j-symbols
          jjjxM = wigner3j(L1,Lx(p),L2,-M1,Mx(p),M2);
          if jjjxM==0, continue; end
          
          jjjxK = wigner3j(L1,Lx(p),L2,-K1,K1-K2,K2);
          if jjjxK==0, continue; end
          
          % Accumulate value
          valp_ = valp_ +  jjjxM * jjjxK * X(p);
        end
        val_ = val_ + (-1)^(K1-M1) * sqrt((2*L1+1)*(2*L2+1)) * valp_;
      end
      
      % Store non-zero value
      if val_==0, continue; end
      idx = idx + 1;
      row(idx) = b1;
      col(idx) = b2;
      values(idx)  = val_;
      if b1~=b2
        idx = idx + 1;
        row(idx) = b2;
        col(idx) = b1;
        values(idx) = val_;
      end
      
    end
  end
  Gamma = sparse(row,col,values,nBasis,nBasis);
    
elseif Method==2
  
  % Method 2: Calculate angular-momentum operator matrices Jp, Jm, Jz, and from
  % them the diffusion operator, for a diagonal diffusion tensor. Calculate
  % potential-dependent term using X coefficients and LMK matrix reps for D^L_MK.
  %-----------------------------------------------------------------------------
  
  % Calculate potential-indepependent part using angular-momentum matrices.
  [Jz,Jp,Jm,J2] = angmomops_LMK(L,M,K);
  Gamma = Rd*(Jp^2+Jm^2) + Rperp*(J2-Jz^2) + Rz*Jz^2;
  
  % Calculate potential-dependent part using D^L_MK matrix representations.
  if usePotential
    DLMK = wignerops_LMK(L,M,K,Lx,Mx,Kx);
    for p = 1:numel(X)
      Gamma = Gamma + X(p)*DLMK{p};
    end
  end

else
  
  error('There is no method %d.',Method);
  
end

return

%===============================================================================
% Calculate matrix representations of angular-momentum operators in LMK basis.
function [Jz,Jp,Jm,J2] = angmomops_LMK(L,M,K)

nBasis = numel(L);

row = zeros(nBasis*3,1);
col = zeros(nBasis*3,1);
val = zeros(nBasis*3,1);
idx = 0;
for b1 = 1:nBasis
  L1 = L(b1);
  M1 = M(b1);
  K1 = K(b1);
  for b2 = 1:nBasis
    L2 = L(b2);
    if L1~=L2, continue; end
    M2 = M(b2);
    if M1~=M2, continue; end
    K2 = K(b2);
    if K1~=K2+1, continue; end
    idx = idx + 1;
    row(idx) = b1;
    col(idx) = b2;
    val(idx) = sqrt(L2*(L2+1)-K2*(K2+1));
  end
end

Jp = sparse(row(1:idx),col(1:idx),val(1:idx),nBasis,nBasis);
Jm = Jp';

didx = 1:nBasis;
J2 = sparse(didx,didx,L.*(L+1));
Jz = sparse(didx,didx,K);

return


%===============================================================================
% Calculate matrix elements of Wigner functions
%
%    <L1,M1,K1|D^L_MK|L2,M2,K2>
%
% D = angmomops(L,M,K,Lx,Mx,Kx)
%
% Inputs:
%   L, K, M:  vectors of L1,M1,K1 and L2,M2,K2
%   Lx, Mx, Kx, vectors of L,M,K for the operator
%
% Outputs:
%   D: cell array of operator matrices, one for each
%      element in Lx,Mx,Kx

function DLMK = wignerops_LMK(L,M,K,Lx,Mx,Kx)

nOps = numel(Lx);
nBasis = numel(L);
for iOp = 1:nOps
  Lx_ = Lx(iOp);
  Mx_ = Mx(iOp);
  Kx_ = Kx(iOp);
  
  % Special case: Lx = 0 (with Mx=Kx=0)
  if Lx_==0
    DLMK{iOp} = speye(nBasis,nBasis);
    continue
  end
  
  D_ = zeros(nBasis,nBasis);
  
  for b1 = 1:nBasis
    L1 = L(b1);
    M1 = M(b1);
    K1 = K(b1);
    for b2 = 1:nBasis
      L2 = L(b2);
      M2 = M(b2);
      K2 = K(b2);
      
      % Screen using selection rules for 3j symbols
      dM = M1-M2;
      dK = K1-K2;
      if Mx_~=dM, continue; end
      if Kx_~=dK, continue; end
      if abs(L1-Lx_)>L2 || L2>L1+Lx_, continue; end
      
      % Calculate 3j symbols, abort as soon as a zero is encountered
      v = wigner3j([L1 Lx_ L2],[-M1 dM M2]);
      if v==0, continue; end
      v = wigner3j([L1 Lx_ L2],[-K1 dK K2])*v;
      if v==0, continue; end
      
      % Evaluate and store non-zero matrix element
      D_(b1,b2) = sqrt((2*L1+1)*(2*L2+1))*(-1)^(K1-M1)*v;
      
    end
  end
  
  DLMK{iOp} = sparse(D_);
  
end
