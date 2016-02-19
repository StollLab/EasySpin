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
%     [1]          simple FID, 3p-ESEEM, echo transient, DEFENCE
%     [1 1]        2p-ESEEM, CP
%     [1 -1]       PEANUT
%     [1 2 1]      2D 3p-ESEEM
%     [1 2]        HYSCORE, DONUT-HYSCORE
%     [1 2 2 1]    2D CP
%     [1 2 -2 1]   2D PEANUT
%     [1 -1 1 -1]
%
%   [1] is the default. For an explanation of
%   the format, see the documentation.
%
%   Mix is a cell array containing the propagators
%   of the mixing sequence(s), for experiments with more
%   than 1 sweep period.
%
%   td is a vector/matrix of the signal with t1 along dim
%   1 and t2 along dim 2.

function Signal = evolve(Sig,Det,Ham,n,dt,IncScheme,Mix)

if (nargin==0), help(mfilename); return; end

if (nargout<1), error('Not enough output arguments.'); end
if (nargout>1), error('Too many output arguments.'); end
if (nargin<5) || (nargin>7), error('Wrong number of input arguments!'); end

if (nargin<6), IncScheme = 1; end
if (nargin<7), Mix = {}; end

if numel(n)~=1 || mod(n,1)~=0 || n<=0
  error('n, the number of points (4th argument), must be a positive integer.');
end

% IncScheme check
%------------------------------------------------------------
if (length(IncScheme)>1) && (nargin<7),
  error('The requested IncScheme requires mixing propagators, but none are provided!');
end
if any((abs(IncScheme)~=1) & (abs(IncScheme)~=2))
  error('IncScheme can contain only 1, -1, 2, and -2.');
end

nEvolutionPeriods = numel(IncScheme);
nDims = max(abs(IncScheme));

% Parameter parsing
%------------------------------------------------------------
if ~iscell(Mix)
  Mix = {Mix};
end
nMix = numel(Mix);

if (nMix~=nEvolutionPeriods-1),
  error('Number of mixing propagators not correct! %d are needed.',nEvolutionPeriods-1);
end
N = size(Sig,1);

if (nDims==1)
  Signal = zeros(n,1);
  if iscell(Ham), Ham = Ham{1}; end
else
  if numel(dt)==1, dt = [dt dt]; end
  if numel(n)==1, n = [n n]; end
  Signal = zeros(n);
end

% Transform all operators in propagator eigenbasis(eigenbases)
%------------------------------------------------------------
if (nDims==1), % 1D case
  % diagonalize propagator
  [Vecs,E] = eig(Ham); % MHz, E doesn't have to be sorted
  E = real(diag(E));
  diagU = exp(-2i*pi*dt*E);
  
  % transform all other matrices to propagator eigenbasis
  Density = Vecs'*Sig*Vecs;
  for iMix = 1:nMix
    Mix{iMix} = Vecs'*Mix{iMix}*Vecs;
  end
  Detector = Vecs'*Det*Vecs;

else % 2D case
  if ~iscell(Ham)
    % Diagonalize propagator
    [Vecs,E] = eig(Ham); % E doesn't have to be sorted
    E = real(diag(E));
    diagUX = exp(-2i*pi*dt(1)*E);
    diagUY = exp(-2i*pi*dt(2)*E);
    
    % Transform all other matrices to propagator eigenbasis
    Density = Vecs'*Sig*Vecs;
    for iMix = 1:nMix
      Mix{iMix} = Vecs'*Mix{iMix}*Vecs;
    end
    Detector = Vecs'*Det*Vecs;
  else
    % Diagonalize propagators
    [V,Ex] = eig(Ham{1}); Vecs{1} = V;
    [V,Ey] = eig(Ham{2}); Vecs{2} = V;
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
end


% Time-domain evolution, IncScheme switchyard
%------------------------------------------------------------
% pre-reshape for trace calculation
Detector = reshape(Detector.',1,N^2);

if isequal(IncScheme,[1]) % IncScheme [1]
  FinalDensity = Density(:);
  UUt = diagU*diagU';
  UUt = UUt(:);
  for ix = 1:n
    Signal(ix) = Detector*FinalDensity;
    FinalDensity = UUt.*FinalDensity;
  end
  
elseif isequal(IncScheme,[1 1]) % IncScheme [1 1]
  % The propagator is diagonal, so we can rewrite
  % U*Mix*U = (diagU*diagU.').*Mix.
  UU = diagU*diagU.';
  % Now we are evolving one dimension. It is not
  % necessary to evolve the initial density matrix.
  % Only the mixing propagator needs to be changed!
  Mix1 = Mix{1};
  for ix = 1:n
    % compute density right before detection
    FinalDensity = Mix1*Density*Mix1';
    % compute trace(Detector*FinalDensity)
    Signal(ix) = Detector*FinalDensity(:);
    Mix1 = UU.*Mix1; % equivalent to U*Mix1*U
  end
  
elseif isequal(IncScheme,[1 -1]) % IncScheme [1 -1]
  MixX = diag(diagU.^n)*Mix{1}; % pre-propagate to end of second period
  UtU = conj(diagU)*diagU.';
  for ix = 1:n
    FinalDensity = MixX*Density*MixX';
    Signal(ix) = Detector*FinalDensity(:);
    MixX = UtU.*MixX;
  end
  
elseif isequal(IncScheme,[1 2]) % IncScheme [1 2]
  UUtX = diagUX*diagUX';
  UUtY = diagUY*diagUY';
  UUtY = reshape(UUtY,N^2,1);
  Mix1 = Mix{1};
  for ix = 1:n(1)
    FinalDensity = reshape(Mix1*Density*Mix1',N^2,1);
    for iy = 1:n(2)
      Signal(ix,iy) = Detector*FinalDensity;
      FinalDensity = UUtY.*FinalDensity;
    end
    Density = UUtX.*Density;
  end
  
elseif isequal(IncScheme,[1 2 1]) % IncScheme [1 2 1]
  Mix1 = Mix{1};
  Mix2 = Mix{2};
  UUX = diagUX*diagUX.';
  UY = diag(diagUY);
  for iy = 1:n(2)
    MixY = Mix2*Mix1;
    MixYt = MixY';
    for ix = 1:n(1)
      FinalDensity = MixY*Density*MixYt;
      Signal(ix,iy) = Detector*FinalDensity(:);
      MixY = UUX.*MixY;
    end
    Mix1 = UY*Mix1;
  end
  
elseif isequal(IncScheme,[1 2 2 1]) % IncScheme [1 2 2 1]
  Mix1 = Mix{1};
  Mix2 = Mix{2};
  Mix3 = Mix{3};
  UUX = diagUX*diagUX.';
  UUY = diagUY*diagUY.';
  for iy = 1:n(2)
    MixY = Mix3*Mix2*Mix1;
    MixYt = MixY';
    for ix = 1:n(1)
      FinalDensity = MixY*Density*MixYt;
      Signal(ix,iy) = Detector*FinalDensity(:);
      MixY = UUX.*MixY;
    end
    Mix2 = UUY.*Mix2;
  end
  
elseif isequal(IncScheme,[1 2 -2 1]) % IncScheme [1 2 -2 1]
  Mix1 = Mix{1};
  Mix2 = diag(diagUY.^n(2))*Mix{2}; % pre-propagate to endpoint of third delay
  Mix3 = Mix{3};
  UUX = diagUX*diagUX.';
  UtUY = conj(diagUY)*diagUY.';
  for iy = 1:n(2)
    MixY = Mix3*Mix2*Mix1;
    MixYt = MixY';
    for ix = 1:n(1)
      FinalDensity = MixY*Density*MixYt;
      Signal(ix,iy) = Detector*FinalDensity(:);
      MixY = UUX.*MixY;
    end
    Mix2 = UtUY.*Mix2;
  end
  
elseif isequal(IncScheme,[1 -1 1 -1]) % IncScheme [1 -1 1 -1]
  Mix1X = diag(diagU.^n)*Mix{1}; % pre-propagate to endpoint of second delay
  Mix2 = Mix{2};
  Mix3X = diag(diagU.^n)*Mix{3}; % pre-propagate to endpoint of fourth delay
  UtU = conj(diagU)*diagU.'; % propagator for Mix1 and Mix3 (add before, remove after)
  for ix = 1:n
    MixX = Mix3X*Mix2*Mix1X;
    FinalDensity = MixX*Density*MixX';
    Signal(ix) = Detector*FinalDensity(:);
    Mix1X = UtU.*Mix1X;
    Mix3X = UtU.*Mix3X;
  end

else
  error('Unsupported incrementation scheme!');
end

return
