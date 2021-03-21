% Calculates the vector representation of sqrt(Peq) in a LMK or LMKjK basis
% under an orientational potential.
%
% Inputs:
%   basis      orientational basis (structure with
%                fields L, M, K, and possibly jK)
%   Potential  orientational potential parameters (structure with fields
%                L, M, K, and lambda; or n x 4 array)
%   Opt        structure with options
%     .useSelectionRules
%     .PeqTolerances
%
% Outputs:
%   sqrtPeq    vector representation of sqrt(Peq)
%   nIntegrals number of evaluated 1D/2D/3D integrals

function [sqrtPeq,nIntegrals] = chili_eqpopvec2(basis,Potential,Opt)

% Parse options
if nargin<3
  Opt = struct;
end
if ~isfield(Opt,'useSelectionRules')
  Opt.useSelectionRules = true;
end
if ~isfield(Opt,'PeqTolerances') || isempty(Opt.PeqTolerances)
  Opt.PeqTolerances = [1e-6 1e-4 1e-4];
end
useSelectionRules = Opt.useSelectionRules;
PeqTolerances = Opt.PeqTolerances;
PeqIntThreshold = PeqTolerances(1);
PeqIntAbsTol = PeqTolerances(2);
PeqIntRelTol = PeqTolerances(3);

% Parse potential
if isstruct(Potential)
  Lp = Potential.L;
  Mp = Potential.M;
  Kp = Potential.K;
  lambda = Potential.lambda;
elseif isnumeric(Potential)
  if size(Potential,2)~=4
    error('Potential must have four columns (L, M, K, lambda)');
  end
  Lp = Potential(:,1);
  Mp = Potential(:,2);
  Kp = Potential(:,3);
  lambda = Potential(:,4);
end

% Process potential
%--------------------------------------------------------------------------
% Assure that potential contains only terms with nonnegative K, and nonnegative M
% for K=0. The other terms needed to render the potential real-valued are
% implicitly supplemented in the subfunction U(a,b,c).
if any(Kp<0)
  error('Only potential terms with nonnegative values of K are allowed. Terms with negative K required to render the potential real-valued are supplemented automatically.');
end
if any(Mp(Kp==0)<0)
  error('For potential terms with K=0, M must be nonnegative. Terms with negative M required to render the potential real-valued are supplemented automatically.');
end
zeroMK = Mp==0 & Kp==0;
if any(~isreal(lambda(zeroMK)))
  error('Potential coefficients for M=K=0 must be real-valued.');
end

% Remove zero entries and constant offset
rmv = lambda==0;
rmv = rmv | (Lp==0 & Mp==0 & Kp==0);
if any(rmv)
  lambda(rmv) = [];
  Lp(rmv) = [];
  Mp(rmv) = [];
  Kp(rmv) = [];
end

% Counter for the number of numerical 1D, 2D, and 3D integrals evaluated.
nIntegrals = [0 0 0];

jKbasis = isfield(basis,'jK') && ~isempty(basis.jK) && any(basis.jK);

% Parse basis
L = basis.L;
M = basis.M;
K = basis.K;
if jKbasis
  jK = basis.jK;
end
nOriBasis = numel(L);

% Handle special case of no potential
if ~any(lambda)
  idx0 = find(L==0 & M==0 & K==0);
  if numel(idx0)~=1
    error('Exactly one orientational basis function with L=M=K=0 is allowed.');
  end
  sqrtPeq = sparse(idx0,1,1,nOriBasis,1);
  return
end

% Abbreviations for 1D, 2D, and 3D integrals
int_b = @(f) integral(f,0,pi,'AbsTol',PeqIntAbsTol,'RelTol',PeqIntRelTol);
int_ab = @(f) integral2(f,0,2*pi,0,pi,'AbsTol',PeqIntAbsTol,'RelTol',PeqIntRelTol);
int_bc = @(f) integral2(f,0,pi,0,2*pi,'AbsTol',PeqIntAbsTol,'RelTol',PeqIntRelTol);
int_abc = @(f) integral3(f,0,2*pi,0,pi,0,2*pi,'AbsTol',PeqIntAbsTol,'RelTol',PeqIntRelTol);

% Calculate partition sum Z for Peq = exp(-U)/Z
zeroMp = all(Mp==0);
zeroKp = all(Kp==0);
if zeroMp && zeroKp
  Z = (2*pi)^2 * int_b(@(b) exp(-U(0,b,0)).*sin(b));
elseif zeroMp && ~zeroKp
  Z = (2*pi) * int_bc(@(b,c) exp(-U(0,b,c)).*sin(b));
elseif ~zeroMp && zeroKp
  Z = (2*pi) * int_ab(@(a,b) exp(-U(a,b,0)).*sin(b));
else
  Z = int_abc(@(a,b,c) exp(-U(a,b,c)).*sin(b));
end
sqrtZ = sqrt(Z);

% Calculate all Wigner coefficients up to Lmax
c = fftso3(@(a,b,c)exp(-U(a,b,c)/2)/sqrtZ,max(L),'DH');

% Rescale
for iL = 1:numel(c)
  c{iL} = c{iL}/sqrt((2*(iL-1)+1)/(8*pi^2));
end

% Copy coefficients into vector
sqrtPeq = zeros(nOriBasis,1);
idxM = M+L+1;
idxL = K+L+1;
for b = 1:nOriBasis
  L_ = L(b);
  sqrtPeq(b) = c{L_+1}(idxM(b),idxL(b));
end

% Remove imaginary part and convert to sparse form
sqrtPeq = real(sqrtPeq);
sqrtPeq(abs(sqrtPeq)<PeqIntThreshold) = 0;
sqrtPeq = sparse(sqrtPeq);

  % General orientational potential function (real-valued)
  % (assumes nonnegative K, and nonnegative M for K=0, and real lambda for
  % M=K=0; other terms are implicitly supplemented)
  function u = U(a,b,c)
    u = 0;
    for p = 1:numel(lambda)
      if lambda(p)==0, continue; end
      if Kp(p)==0 && Mp(p)==0
        u = u - wignerd([Lp(p) Mp(p) Kp(p)],b) * real(lambda(p)) - 0*a;
      else
        u = u - 2*real(wignerd([Lp(p) Mp(p) Kp(p)],a,b,c) * lambda(p));
      end
    end
  end

end
