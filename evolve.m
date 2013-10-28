% evolve  Time domain evolution of density matrix
%
%   td = evolve(Sig,Det,Ham,n,dt);
%   td = evolve(Sig,Det,Ham,n,dt,IncScheme);
%   td = evolve(Sig,Det,Ham,n,dt,IncScheme,Mix);
%
%   Evolves the density matrix Sig  under the
%   Hamiltonian Ham with time step dt n-1 times
%   and detects using Det after each step. Hermitian
%   input matrices are assumed. td(1) is the value
%   obtained by detecting Sig without evolution.
%
%   IncScheme determines the incrementation scheme
%   and can be one of the following (up to four
%   sweep periods, up to two dimensions)
%
%     [1]         simple FID, 3p, DEFENCE
%     [1 1]       2p, CP
%     [1 -1]      PEANUT
%     [1 2 1]     2D 3p
%     [1 2]       HYSCORE, DONUT-HYSCORE
%     [1 2 2 1]   2D CP
%     [1 2 -2 1]  2D PEANUT
%
%   [1] is the default. For an explanation of
%   the format, see the documentation.
%
%   Mix is a cell array containing the propagators
%   of the mixing sequence(s), for experiments with more than
%   1 sweep period.
%
%   td is a vector/matrix with t1 along dim
%   1 and t2 along dim 2.

function Signal = evolve(Sig,Det,Ham,n,dt,IncScheme,Mix)

if (nargin==0), help(mfilename); return; end

if (nargout<1), error('Not enough output arguments.'); end
if (nargout>1), error('Too many output arguments.'); end
if (nargin<5) || (nargin>7), error('Wrong number of input arguments!'); end

if (nargin<6), IncScheme = 1; end

if numel(n)~=1 || mod(n,1)~=0 || n<=0
  error('n, the number of points (4th argument), must be a positive integer.');
end

% IncScheme check
%------------------------------------------------------------
ID = 0;

switch numel(IncScheme)
case 1
   if IncScheme==1, ID = 1; end
case 2
  if all(IncScheme==[1 1]), ID = 2;
  elseif all(IncScheme==[1 -1]), ID = 3;
  elseif all(IncScheme==[1 2]), ID = 4;
  end
case 3
  if all(IncScheme==[1 2 1]), ID = 5; end
case 4
  if all(IncScheme==[1 2 2 1]), ID = 6;
  elseif all(IncScheme==[1 2 -2 1]), ID = 7;
  end
end

if ~ID,
  error('Unsuppported incrementation scheme!');
end
if (length(IncScheme)>1) && (nargin<7),
  error('The requested IncScheme requires mixing propagators, but the 7th input argument is missing!');
end

% Parameter parsing
%------------------------------------------------------------
if (nargin<7)
  Mix = [];
end
if (~iscell(Mix))
  if isempty(Mix)
    Mix = {};
  else
    if (ndims(Mix)==2)
      Mix = {Mix};
    else
      for m=1:size(Mix,3)
        M{m} = Mix(:,:,m);
      end
      Mix = M;
    end
  end
end
nMix = numel(Mix);

if (nMix~=length(IncScheme)-1),
  error('Number of mixing propagators not correct!');
end
nDims = 1 + (ID>3);
N = size(Sig,1);
NN = N^2;

% Transform all operators in propagator eigenbasis(eigenbases)
%------------------------------------------------------------
if nDims==1, % 1D case
  % diagonalize propagator
  [Vecs,E] = eig(Ham); % MHz, E doesn't have to be sorted
  diagU = exp(-2i*pi*dt*diag(real(E)));
  
  % transform all other matrices
  Density = Vecs'*Sig*Vecs;
  for iMix = 1:nMix
    Mix{iMix} = Vecs'*Mix{iMix}*Vecs;
  end
  Detector = Vecs'*Det*Vecs;
  Signal = zeros(n,1);

else %2D case
  if numel(dt)==1, dt = [dt dt]; end
  if numel(n)==1, n = [n n]; end
  if size(Ham,3)==1,
    % Diagonalize propagator
    [Vecs,E] = eig(Ham); % E doesn't have to be sorted
    E = real(diag(E));
    diagUX = exp(-2i*pi*dt(1)*E);
    diagUY = exp(-2i*pi*dt(2)*E);
    
    % Transform all other matrices
    Density = Vecs'*Sig*Vecs;
    for iMix = 1:nMix
      Mix{iMix} = Vecs'*Mix{iMix}*Vecs;
    end
    Detector = Vecs'*Det*Vecs;
  else
    % Diagonalize propagators
    [V,Ex] = eig(Ham(:,:,1)); Vecs{1} = V;
    [V,Ey] = eig(Ham(:,:,2)); Vecs{2} = V;
    diagUX = exp((-2i*pi*dt(1))*real(diag(Ex)));
    diagUY = exp((-2i*pi*dt(2))*real(diag(Ey)));

    % Transform all other matrices
    d = abs(IncScheme);
    Density = Vecs{d(1)}'*Sig*Vecs{d(1)};
    for iMix = 1:nMix
      Mix{iMix} = Vecs{d(iMix+1)}'*Mix{iMix}*Vecs{d(iMix)};
    end
    Detector = Vecs{d(end)}'*Det*Vecs{d(end)};
  end
  Signal = zeros(n);
end


% Time-domain evolution, IncScheme switchyard
%------------------------------------------------------------
% pre-reshape for trace calculation
Detector = reshape(Detector.',1,NN);

switch ID

case 1 % IncScheme [1]
  FinalDensity = Density(:);
  UUt = diagU*diagU';
  UUt = UUt(:);
  for k = 1:n
    Signal(k) = Detector*FinalDensity;
    FinalDensity = UUt.*FinalDensity;
  end
  
case 2 % IncScheme [1 1]
  % The propagator is diagonal, so we can rewrite
  % U*Mix*U = (diagU*diagU.').*Mix.
  UU = diagU*diagU.';
  % Now we are evolving one dimension. It is not
  % necessary to evolve the initial density matrix,
  % since we add a new U to Mix both from the left
  % and from the right!!
  Mix1 = Mix{1};
  for k = 1:n
    % compute density right before detection
    FinalDensity = Mix1*Density*Mix1';
    % compute trace(Detector*FinalDensity)
    Signal(k) = Detector*FinalDensity(:);
    Mix1 = UU.*Mix1; % equivalent to U*Mix1*U
  end
  
case 3  % IncScheme [1 -1]
  MixX = diag(diagU.^n)*Mix{1};
  MixXt = MixX';
  UtU = conj(diagU)*diagU.';
  for x = 1:n
    FinalDensity = MixX*Density*MixXt;
    Signal(x) = Detector*FinalDensity(:);
    MixX = UtU.*MixX;
  end
  
case 4 % IncScheme [1 2]
  UUtX = diagUX*diagUX';
  UUtY = diagUY*diagUY';
  UUtY = reshape(UUtY,NN,1);
  Mix1 = Mix{1};
  for x = 1:n(1)
    FinalDensity = reshape(Mix1*Density*Mix1',NN,1);
    for y = 1:n(2)
      Signal(x,y) = Detector*FinalDensity;
      FinalDensity = UUtY.*FinalDensity;
    end
    Density = UUtX.*Density;
  end
  
case 5 % IncScheme [1 2 1]
  Mix1 = Mix{1};
  Mix2 = Mix{2};
  UUX = diagUX*diagUX.';
  UY = diag(diagUY);
  for y = 1:n(2)
    MixY = Mix2*Mix1;
    MixYt = MixY';
    for x = 1:n(1)
      FinalDensity = MixY*Density*MixYt;
      Signal(x,y) = Detector*FinalDensity(:);
      MixY = UUX.*MixY;
    end
    Mix1 = UY*Mix1;
  end
  
case 6 % IncScheme [1 2 2 1]
  Mix1 = Mix{1};
  Mix2 = Mix{2};
  Mix3 = Mix{3};
  UUX = diagUX*diagUX.';
  UUY = diagUY*diagUY.';
  for y = 1:n(2)
    MixY = Mix3*Mix2*Mix1;
    MixYt = MixY';
    for x = 1:n(1)
      FinalDensity = MixY*Density*MixYt;
      Signal(x,y) = Detector*FinalDensity(:);
      MixY = UUX.*MixY;
    end
    Mix2 = UUY.*Mix2;
  end
  
case 7 % IncScheme [1 2 -2 1]
  Mix1 = Mix{1};
  Mix2 = diag(diagUY.^n(2))*Mix{2};
  Mix3 = Mix{3};
  UUX = diagUX*diagUX.';
  UtUY = conj(diagUY)*diagUY.';
  for y = 1:n(2)
    MixY = Mix3*Mix2*Mix1;
    MixYt = MixY';
    for x = 1:n(1)
      FinalDensity = MixY*Density*MixYt;
      Signal(x,y) = Detector*FinalDensity(:);
      MixY = UUX.*MixY;
    end
    Mix2 = UtUY.*Mix2;
  end
end

return
