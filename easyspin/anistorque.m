%  anistorque  Calculate the torque due to an anisotripic orienting potential.
%
%   torque = anistorque(LMK, Coefs, q);
%
%   Input:
%     LMK            numeric, size = (nCoefs,3)
%                    integers corresponding to the quantum numbers L, M, 
%                    and K
%
%     Coefs          numeric, size = (nCoefs,2)
%                    real coefficients cLMK+ and cLMK- for orienting 
%                    potentials
%
%     q              numeric, size = (4,nTraj)
%                    quaternions
%
%   Output:
%     anistorque     numeric, size = (3,1,nTraj)

%   References
%   ----------
%   [1] Sezer, et al., J.Chem.Phys. 128, 165106 (2008), doi: 10.1063/1.2908075

function torque = anistorque(LMK, Coefs, q)
%% Preprocessing
shapeq = num2cell(size(q));
Index = cell(1, ndims(q));
Index(:) = {':'};

if ~isnumeric(q) || size(q,1)~=4
  error('q must be an array of size 4x....')
end

if ~ismatrix(LMK) || size(LMK,2)~=3
  error('LMK must be an array of shape Nx3.')
end

if ~ismatrix(Coefs) || size(Coefs,2)~=2
  error('Coefs must be an array of shape Nx2.')
end

nCoefs = size(Coefs,1);

lambda = Coefs(:,1) + 1i*Coefs(:,2);

Lvals = LMK(:,1);
Mvals = LMK(:,2);
Kvals = LMK(:,3);

if any(Lvals(:)<1)
  error('All values of L must be greater than or equal to one.')
end

if any(Mvals(:)<-Lvals)
  error('All values of M must be greater than or equal to -L.')
end

if any(Kvals(:)>Lvals)
  error('All values of K must be less than or equal to L.')
end

iJvecu = zeros(3,size(q,2),nCoefs);

%%
% Not sure if there is a good way to vectorize this
for n=1:nCoefs
  L = Lvals(n);
  M = Mvals(n);
  K = Kvals(n);
  iJvecu(Index{:},n) = real(iJu(L,M,K,lambda(n),q));
%  Why doesn't this work?! Accumulation of rounding errors?
%   iJvecu(Index{:},n) =               0.5*(iJu(L,M,K,lambda(n),q) ...
%                               + (-1)^(M-K)*iJu(L,-M,-K,conj(lambda(n)),q));
end

if nCoefs > 1
  torque = -sum(iJvecu,ndims(q)+1);
else
  torque = -iJvecu;
end

  assert(all(imag(torque(:))<1e-14), 'Torque is not real.')

%% Helper functions
function  out = iJu(L, M, K, cLK, q)
% Calculate the result of i\vo{J}u for a given potential coefficient
% cjm corresponds to c_{m}^{j} in Eq. 58 in [1]

if K==L
  % K cannot be greater than L, so D_{0,K+1}^L = 0
  iJxu = 1i/2*cLK*Cpm(L,K,'-')*wignerdquat(L,M,K-1,q);
  iJyu = -1/2*cLK*Cpm(L,K,'-')*wignerdquat(L,M,K-1,q);
elseif K==-L
  % K cannot be less than -L, so D_{0,K-1}^L = 0
  iJxu = 1i/2*cLK*Cpm(L,K,'+')*wignerdquat(L,M,K+1,q);
  iJyu =  1/2*cLK*Cpm(L,K,'+')*wignerdquat(L,M,K+1,q);
else
  iJxu = 1i/2*cLK*(Cpm(L,K,'+')*wignerdquat(L,M,K+1,q) ...
                  +Cpm(L,K,'-')*wignerdquat(L,M,K-1,q));
  iJyu =  1/2*cLK*(Cpm(L,K,'+')*wignerdquat(L,M,K+1,q) ...
                  -Cpm(L,K,'-')*wignerdquat(L,M,K-1,q));
end

iJzu = 1i*cLK*K*wignerdquat(L,M,K,q);    

out = [iJxu; iJyu; iJzu];

end


function C = Cpm(L, K, pm)
% Eq. 60 in [1]

if     pm=='+', C = sqrt(L*(L+1)-K*(K+1));
elseif pm=='-', C = sqrt(L*(L+1)-K*(K-1)); end

end

end
