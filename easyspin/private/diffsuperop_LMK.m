% This function computes the diffusion superoperator matrix in the LMK basis,
% for a diagonal or non-diagonal diffusion tensor. It is used by chili().
%
% It does not support an ordering potential.
% It does not require any particular order of orientational basis functions.
%
% Syntax:
%   Gamma = diffsuperop_LMK(basis,R,XLK,Method)
%
% Input:
%   basis    structure containing the lists for the quantum numbers for
%            the orientational basis in basis.L, basis.M, and basis.K
%   R        array with 3 principal values of diffusion tensor (s^-1),
%            or single value if diffusion is isotropic, or full 3x3 matrix
%   XLK      (optional)
%            ordering potential coefficients, as returned by
%            chili_xlk (assumed all zero if not given)
%   Method   (optional)
%            1: method for diagonal diffusion tensors
%            2: method for diagonal diffusion tensors (a little slower than 1)
%            3: method for 3x3 diffusion tensors (needed for non-diagonal R)
%
% Output:
%   Gamma    diffusion superoperator matrix in the given LMK basis (s^-1),
%            sparse

function Gamma = diffsuperop_LMK(basis,R,XLK,Method)

if isfield(basis,'jK') && ~isempty(basis.jK)
  error('This function expects an LMK basis, without jK.');
end

if nargin<4, Method = []; end

usePotential = nargin==3 && ~isempty(XLK) && any(XLK(:)~=0);
if usePotential
  Lxmax = size(XLK,1)-1;
end

fullDiffusionTensor = false;
switch numel(R)
  case 1 % isotropic value
    if isempty(Method), Method = 1; end
    Rd = 0;
    Rperp = R;
    Rz = R;
  case 3 % 3 principal values
    if isempty(Method), Method = 1; end
    Rd = (R(1)-R(2))/4;
    Rperp = (R(1)+R(2))/2;
    Rz = R(3);
  case 9 % full tensor
    fullDiffusionTensor = true;
    if isempty(Method), Method = 3; end
  otherwise
    error('Wrong number of elements in diffusion tensor.');
end

% Calculate diffusion operator
%-------------------------------------------------------------------------------
L = basis.L;
M = basis.M;
K = basis.K;
nBasis = numel(L);

% Treat the cases of isotropic and axial diffusion tensors in the absence of
% an ordering potential. In these cases, the diffusion operator matrix is
% diagonal.
if ~fullDiffusionTensor && Rd==0 && ~usePotential
  diagonal = Rperp*(L.*(L+1)-K.^2) + Rz*K.^2;
  Gamma = spdiags(diagonal,0,nBasis,nBasis);
  return
end

if Method==1
  % Method 1: Calculate the diffusion operator directly (fast).
  %-----------------------------------------------------------------------------
  % Only diagonal and upper triangular part of matrix are explicitly calulcated,
  % since Gamma is symmetric
  if fullDiffusionTensor
    error('Cannot use full difussion tensor with this method.');
  end
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
        for Lx = abs(L1-L2):min(Lxmax,L1+L2)
          idx_xL = Lx+1;
          if ~any(XLK(idx_xL,:)), continue; end
          if abs(K1-K2)>Lx, continue; end
          
          % Get potential coefficient
          xlk_ = XLK(idx_xL,(K1-K2)+Lx+1);
          if xlk_==0, continue; end
          
          % Calculate M- and K-dependent 3j-symbols
          if abs(M1)>L1 || abs(M1)>L2, continue; end
          jjjxM = wigner3j(L1,Lx,L2,M1,0,-M1);
          if abs(K1)>L1 || abs(K2-K1)>Lx || abs(K2)>L2, continue; end
          jjjxK = wigner3j(L1,Lx,L2,K1,K2-K1,-K2);
          
          % Accumulate value
          valp_ = valp_ +  jjjxM * jjjxK * xlk_;
        end
        if valp_~=0
          val_ = val_ + (-1)^(K1-M1) * sqrt((2*L1+1)*(2*L2+1)) * valp_;
        end
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
  % them the diffusion operator, for a diagonal diffusion tensor.
  %-----------------------------------------------------------------------------
  if fullDiffusionTensor
    error('Cannot use full difussion tensor with this method.');
  end
  if usePotential
    error('Cannot use orienting potential with this method.');
  end
  
  [Jp,Jm,Jz,J2] = angmomops(L,M,K);
  
  Gamma = Rperp*(J2-Jz^2) + Rz*Jz^2 + Rd*(Jp^2+Jm^2);
  
elseif Method==3
  
  % Method 3: Calculate angular-momentum operator matrices Jx, Jy, Jz, and from
  % them the diffusion operator, for a diagonal or non-diagonal diffusion tensor.
  %-----------------------------------------------------------------------------
  if usePotential
    error('Cannot use orienting potential with this method.');
  end
  
  [Jp,Jm,Jz] = angmomops(L,M,K);
  Jx = (Jp+Jm)/2;
  Jy = (Jp-Jm)/2i;
  
  if numel(R)==3
    Gamma = R(1)*Jx^2 + R(2)*Jy^2 + R(3)*Jz^2;
  else
    Gamma = sparse(nBasis,nBasis);
    J = {Jx,Jy,Jz};
    for i = 1:3
      for j = 1:3
        Gamma = Gamma + R(i,j)*J{i}*J{j};
      end
    end
  end
  
else
  
  error('Only Methods 1, 2, and 3 are implemented.');
  
end

return

%===============================================================================
% Calculate matrix representations of common angular-momentum operators in
% LMK basis.
function [Jp,Jm,Jz,J2] = angmomops(L,M,K)

nBasis = numel(L);

Jp = zeros(nBasis,nBasis);
Jm = zeros(nBasis,nBasis);
Jz = zeros(nBasis,nBasis);
J2 = zeros(nBasis,nBasis);
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
    if K1==K2
      Jz(b1,b2) = K2;
      J2(b1,b2) = L2*(L2+1);
    elseif K1==K2+1
      Jp(b1,b2) = sqrt(L2*(L2+1)-K2*(K2+1));
    elseif K1==K2-1
      Jm(b1,b2) = sqrt(L2*(L2+1)-K2*(K2-1));
    end
  end
end
Jp = sparse(Jp);
Jm = sparse(Jm);
Jz = sparse(Jz);
J2 = sparse(J2);

return
