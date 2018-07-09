function stvec = chili_startingvector(Basis,Potential)

if any(Potential.lambda) && ~Potential.oldstyle
  error('This functions works only for potentials with M=0, L=2,4, and K=0,2.');
end

lambda = Potential.lambda;
isPotential = any(lambda);
maxKp = max(Potential.K);

isNuc1 = isfield(Basis,'pI1');
isNuc2 = isfield(Basis,'pI2');

nBasis = numel(Basis.L);
stvec = zeros(nBasis,1);

absTol = 1e-8; % for numerical integration
useOldIntegrator = verLessThan('matlab','7.14'); % us quadl() instead of integral()

cacheval = NaN;
cacheLK = [NaN NaN];

for b = 1:nBasis
  L = Basis.L(b);
  K = Basis.K(b);
  
  if mod(L,2) || mod(K,2), continue; end
  if Basis.M(b)~=0, continue; end
  if Basis.jK(b)~=1, continue; end
  if abs(Basis.pS(b))~=1, continue; end
  if isNuc1 && Basis.pI1(b)~=0, continue; end
  if isNuc2 && Basis.pI2(b)~=0, continue; end
  
  if ~isPotential
    % Zero potential: non-zero value only for L==K==0
    if L~=0 || K~=0, continue; end
    stvec(b) = 1;
  else
    % Zero for axial potential (max Kp == 0) and K>0
    if maxKp==0 && K~=0, continue; end
    % Numerical integration if value is different from previous one
    if all([L K]==cacheLK)
      value = cacheval;
    else
      fun  = @(z)orifun(z,L,K,lambda);
      if useOldIntegrator
        value = quadl(fun,0,1,absTol); %#ok<DQUADL>
      else
        value = integral(fun,0,1,'AbsTol',absTol);
      end
      value = value * sqrt((2*L+1)/prod(L-K+1:L+K));
      if K==0
        value = value/sqrt(2);
      end
      cacheval = value;
      cacheLK = [L K];
    end
    stvec(b) = value;
  end
  
end

stvec = stvec/norm(stvec);
stvec = sparse(stvec);

return

%===============================================================================
% Integrand of the orientational integral in the starting vector.
% (see Schneider 1989, p.19)
function val = orifun(z,L,K,lambda)

% based on the integral (valid for even K)
%   \int_0^{2\pi} cos(K g) exp(B cos(2 g)) d g = 2 \pi besseli(K/2,B)

% Potential terms with K == 0
A = 0;
if lambda(1), A = A + lambda(1)/2*plegendre(2,0,z); end
if lambda(3), A = A + lambda(3)/2*plegendre(4,0,z); end

% Potential terms with K == 2
B = 0;
if lambda(2), B = B + lambda(2)/(2*sqrt( 6))*plegendre(2,2,z); end
if lambda(4), B = B + lambda(4)/(6*sqrt(10))*plegendre(4,2,z); end

val = plegendre(L,K,z).*exp(A).*besseli(K/2,B)*(2*pi);

return
