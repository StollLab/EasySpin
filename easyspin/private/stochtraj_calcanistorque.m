%  stochtraj_calcanistorque  Calculate the torque due to an anisotripic orienting potential.
%
%   torque = stochtraj_calanistorque(LMK, lambda, q);
%
%   Input:
%     LMK            numeric, size = (nCoefs,3)
%                    integers corresponding to the quantum numbers L, M, 
%                    and K
%
%     lambda         numeric, size = (nCoefs,1)
%                    real coefficients cLMK+ and cLMK- for orienting 
%                    potentials
%
%     q              numeric, size = (4,nTraj)
%                    quaternions
%
%
%   Output:
%     anistorque     numeric, size = (3,1,nTraj)

%   References
%   ----------
%   [1] Sezer, et al., J.Chem.Phys. 128, 165106 (2008), doi: 10.1063/1.2908075

function torque = stochtraj_calcanistorque(LMK, lambda, q)
% Preprocessing
% -------------------------------------------------------------------------
shapeq = num2cell(size(q));
Index = cell(1, ndims(q));
Index(:) = {':'};

if ~isnumeric(q) || size(q,1)~=4
  error('q must be an array of size 4x....')
end

if ~ismatrix(LMK) || size(LMK,2)~=3
  error('LMK must be an array of shape Nx3.')
end

if ~ismatrix(lambda) || size(lambda,2)~=1
  error('lambda must be an array of shape Nx1.')
end

nCoefs = size(lambda,1);

Lvals = LMK(:,1);
Mvals = LMK(:,2);
Kvals = LMK(:,3);

iJvecu = zeros(3, size(q,2), nCoefs);

% Not sure if there is a good way to vectorize this
for n = 1:nCoefs
  L = Lvals(n);
  M = Mvals(n);
  K = Kvals(n);
  
  if M==0 && K==0
    iJvecu(Index{:},n) = real(iJu(L,M,K,lambda(n),q));
  else
    iJvecu(Index{:},n) = 2*real(iJu(L, M, K, lambda(n), q));
  end
end

if nCoefs > 1
  torque = -sum(iJvecu,ndims(q)+1);
else
  torque = -iJvecu;
end

if all(imag(torque(:))<1e-11)
  torque = real(torque);
else
  error('Finite imaginary part of torque detected.')
end

end
 

% Helper functions
% -------------------------------------------------------------------------
function  out = iJu(L, M, K, lam, q)
% Calculate the result of i\vo{J}u for a given potential coefficient
% cjm corresponds to c_{m}^{j} in Eq. 58 in [1]

if K==L
  % K cannot be greater than L, so D_{0,K+1}^L = 0
  iJxu = 1i/2*lam*Cpm(L,K,'-')*stochtraj_wignerdquat(L,M,K-1,q);
  iJyu = -1/2*lam*Cpm(L,K,'-')*stochtraj_wignerdquat(L,M,K-1,q);
elseif K==-L
  % K cannot be less than -L, so D_{0,K-1}^L = 0
  iJxu = 1i/2*lam*Cpm(L,K,'+')*stochtraj_wignerdquat(L,M,K+1,q);
  iJyu =  1/2*lam*Cpm(L,K,'+')*stochtraj_wignerdquat(L,M,K+1,q);
else
  iJxu = 1i/2*lam*(Cpm(L,K,'+')*stochtraj_wignerdquat(L,M,K+1,q) ...
                  +Cpm(L,K,'-')*stochtraj_wignerdquat(L,M,K-1,q));
  iJyu =  1/2*lam*(Cpm(L,K,'+')*stochtraj_wignerdquat(L,M,K+1,q) ...
                  -Cpm(L,K,'-')*stochtraj_wignerdquat(L,M,K-1,q));
end

iJzu = 1i*lam*K*stochtraj_wignerdquat(L,M,K,q);    

out = [iJxu; iJyu; iJzu];

end


function C = Cpm(L, K, pm)
% Eq. 60 in [1]

if     pm=='+', C = sqrt(L*(L+1)-K*(K+1));
elseif pm=='-', C = sqrt(L*(L+1)-K*(K-1)); end

end