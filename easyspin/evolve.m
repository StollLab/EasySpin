% evolve  Time-domain evolution of density matrix
%
%   td = evolve(Sig,Det,Ham,n,dt);
%   td = evolve(Sig,Det,Ham,n,dt,IncScheme);
%   td = evolve(Sig,Det,Ham,n,dt,IncScheme,Mix);
%
%   Evolves the density matrix Sig under the Hamiltonian Ham with time
%   step dt n-1 times and detects using Det after each step. Hermitian
%   input matrices are assumed. Ham is assumed to be in units of MHz.
%
%   n gives the number of points along each dimension (a number for 1D
%   experiments, a 2-vector for 2D experiments). dt contains the time increment
%   along each dimension, in microseconds.
%
%   IncScheme determines the incrementation scheme and can be one of the
%   following (up to five incrementation periods, up to two dimensions)
%
%     1D experiments:
%       [1]           simple FID, 3p-ESEEM, echo transient, DEFENCE
%       [1 1]         2p-ESEEM, CP, 3p- and 4p-RIDME
%       [1 -1]        3p-DEER, 4p-DEER, PEANUT, 5p-RIDME
%       [1 -1 1 -1]
%       [1 1 -1 -1]   SIFTER
%       [1 -1 -1 1]   7p-DEER
%
%     2D experiments:
%       [1 2]         3p-ESEEM echo transient, HYSCORE, DONUT-HYSCORE
%       [1 1 2]       2p-ESEEM etc. with echo transient
%       [1 -1 2]      3p-DEER, 4p-DEER etc. with echo transient
%       [1 2 1]       2D 3p-ESEEM
%       [1 1 2 2]     2D refocused 2p-echo
%       [1 2 2 1]     2D CP
%       [1 2 -2 1]    2D PEANUT
%       [1 1 -1 -1 2] SIFTER with echo transient
%       [1 -1 -1 1 2] 7p-DEER with echo transient
%
%   [1] is the default. 1 and -1 indicate incrementation and decrementation
%   along the first dimensions, and 2 and -2 are analogous for the second
%   dimensions. For more details, see the documentation.
%
%   Mix is a cell array containing the propagators of the mixing
%   block(s). It is required for experiments with more than 1 sweep period.
%
%   td is the output signal, a vector for 1D experiments and a matrix for
%   2D experiments. td(1) is the value obtained by detecting Sig without evolution.

function Signal = evolve(Sig,Det,Ham,n,dt,IncScheme,Mix)

if nargin==0, help(mfilename); return; end

if nargout<1, error('Not enough output arguments.'); end
if nargout>1, error('Too many output arguments.'); end
if nargin<5 || nargin>7, error('Wrong number of input arguments!'); end

if nargin<6, IncScheme = 1; end
if nargin<7, Mix = {}; end

if any(mod(n,1)~=0) || any(n<=0)
  error('n, the number of points (4th argument), must be a positive integer.');
end

% IncScheme check
%------------------------------------------------------------
if length(IncScheme)>1 && nargin<7
  error('The requested IncScheme requires mixing propagators, but none are provided!');
end
if ~isnumeric(IncScheme) || ~isvector(IncScheme) || ~isreal(IncScheme)
  error('IncScheme (6th input) must be a vector containing 1, -1, 2, and -2.');
end
if any((abs(IncScheme)~=1) & (abs(IncScheme)~=2))
  error('IncScheme (6th input) can contain only 1, -1, 2, and -2.');
end

nEvolutionPeriods = numel(IncScheme);
nDimensions = max(abs(IncScheme));

% Parameter parsing
%------------------------------------------------------------
if ~iscell(Det)
  Det = {Det};
end
nDetectors = numel(Det);

if ~iscell(Mix)
  Mix = {Mix};
end
nMixingBlocks = numel(Mix);

if nMixingBlocks~=nEvolutionPeriods-1
  error('Number of mixing propagators not correct! %d are needed.',nEvolutionPeriods-1);
end
N = size(Sig,1);

if nDimensions==1
  for iDet = 1:nDetectors
    Signal{iDet} = zeros(n,1);
  end
  if iscell(Ham), Ham = Ham{1}; end
else
  if numel(dt)==1, dt = [dt dt]; end
  if numel(n)==1, n = [n n]; end
  for iDet = 1:nDetectors
    Signal{iDet} = zeros(n);
  end
end
if nDetectors==1
  Signal = Signal{1};
end

% Transform all operators to Hamiltonian eigenbasis (eigenbases)
%---------------------------------------------------------------
if ~iscell(Ham)
  
  if nnz(Ham)==nnz(diag(Ham)) % Check if Hamiltonian is already diagonal
    E = diag(Ham);
    Density = Sig;
    Detector = Det;
  else
    % Diagonalize Hamiltonian
    [Vecs,E] = eig(Ham); % MHz, E doesn't have to be sorted
    E = real(diag(E));
    % Transform all other matrices to Hamiltonian eigenbasis
    Density = Vecs'*Sig*Vecs;
    for iMix = 1:nMixingBlocks
      Mix{iMix} = Vecs'*Mix{iMix}*Vecs;
    end
    for iDet = 1:nDetectors
      Detector{iDet} = Vecs'*Det{iDet}*Vecs;
    end
  end
  
  % Define free evolution propagators
  if nDimensions==1
    diagU = exp(-2i*pi*dt*E);
  else
    diagUX = exp(-2i*pi*dt(1)*E);
    diagUY = exp(-2i*pi*dt(2)*E);
  end
  
else
  
  % Check if Hamiltonians are already diagonal
  if (nnz(Ham{1})==nnz(diag(Ham{1})) && nnz(Ham{1})==nnz(diag(Ham{1})))
    Ex = Ham{1};
    Ey = Ham{2};
    Density = Sig;
    Detector = Det;
  else
    
    % Diagonalize Hamiltonians
    [Vecs{1},Ex] = eig(Ham{1});
    [Vecs{2},Ey] = eig(Ham{2});
    
    % Transform all other matrices to Hamiltonian eigenbasis
    d = abs(IncScheme);
    Density = Vecs{d(1)}'*Sig*Vecs{d(1)};
    for iMix = 1:nMixingBlocks
      Mix{iMix} = Vecs{d(iMix+1)}'*Mix{iMix}*Vecs{d(iMix)};
    end
    for iDet = 1:nDetectors
      Detector{iDet} = Vecs{d(end)}'*Det{iDet}*Vecs{d(end)};
    end
    
  end
  
  % Define free evolution propagators
  diagUX = exp((-2i*pi*dt(1))*real(diag(Ex)));
  diagUY = exp((-2i*pi*dt(2))*real(diag(Ey)));
  
end

% Time-domain evolution, IncScheme switchyard
%------------------------------------------------------------
% The following implementations in the propagator eigenframes
% (giving diagonal propagators) make use of the following simplifications
% of the matrix multiplications associated with the propagations:
%
%   U*Density*U'   = (diagU*diagU').*Density
%   U*Propagator*U = (diagU*diagU.').*Propagator
%   (U^-1)*Propagator*U = U'*Propagator*U = (conj(diagU)*diagU.').*Propagator
%   U*Propagator*(U^-1) = U*Propagator*U' = (diagU*diagU').*Propagator
%
% where U are diagonal matrices and diagU are vectors of eigenvalues.

% Pre-reshape for trace calculation
for iDet = 1:nDetectors
  Detector{iDet} = reshape(Detector{iDet}.',1,N^2);
end
if nDetectors==1
  Detector = Detector{1};
end

if isequal(IncScheme,1)
  
  FinalDensity = Density(:);
  U_ = diagU*diagU';
  U_ = U_(:);
  for ix = 1:n
    % Compute trace(Detector*FinalDensity)
    if nDetectors==1
      Signal(ix) = Detector*FinalDensity;
    else
      for iDet = 1:nDetectors
        Signal{iDet}(ix) = Detector{iDet}*FinalDensity;
      end
    end
    FinalDensity = U_.*FinalDensity; % equivalent to U*FinalDensity*U'
  end
  
elseif isequal(IncScheme,[1 1])
  
  UU_ = diagU*diagU.';
  % It is not necessary to evolve the initial density matrix.
  % Only the mixing propagator needs to be evolved.
  Mix1 = Mix{1};
  for ix = 1:n
    % Compute density right before detection
    FinalDensity = Mix1*Density*Mix1';
    % Compute trace(Detector*FinalDensity)
    if nDetectors==1
      Signal(ix) = Detector*FinalDensity(:);
    else
      for iDet = 1:nDetectors
        Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
      end
    end
    Mix1 = UU_.*Mix1; % equivalent to U*Mix1*U
  end
  
elseif isequal(IncScheme,[1 -1])
  
  % Pre-propagate mixing propagator to end of second period (= start of
  % experiment)
  MixX = diag(diagU.^n)*Mix{1};
  UtU_ = conj(diagU)*diagU.';
  for ix = 1:n
    FinalDensity = MixX*Density*MixX';
    if nDetectors==1
      Signal(ix) = Detector*FinalDensity(:);
    else
      for iDet = 1:nDetectors
        Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
      end
    end
    MixX = UtU_.*MixX; % equivalent to U^-1*MixX*U
  end
  
elseif isequal(IncScheme,[1 2])
  
  UX_ = diagUX*diagUX';
  UY_ = diagUY*diagUY';
  UY_ = reshape(UY_,N^2,1);
  Mix1 = Mix{1};
  for ix = 1:n(1)
    FinalDensity = reshape(Mix1*Density*Mix1',N^2,1);
    for iy = 1:n(2)
      if nDetectors==1
        Signal(ix,iy) = Detector*FinalDensity;
      else
        for iDet = 1:nDetectors
          Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity;
        end
      end
      FinalDensity = UY_.*FinalDensity; % equivalent to UY*Densty*UY';
    end
    Density = UX_.*Density; % equivalent to UX*Densty*UX';
  end
  
elseif isequal(IncScheme,[1 1 2])
  
  UUX_ = diagUX*diagUX.';
  UY_ = diagUY*diagUY';
  Mix1 = Mix{1};
  Mix2 = Mix{2};
  for ix = 1:n(1)
    M = Mix2*Mix1;
    FinalDensity = M*Density*M';
    for iy = 1:n(2)
      if nDetectors==1
        Signal(ix,iy) = Detector*FinalDensity(:);
      else
        for iDet = 1:nDetectors
          Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
        end
      end
      FinalDensity = UY_.*FinalDensity; % equivalent to UY*FinalDensity*UY'
    end
    Mix1 = UUX_.*Mix1; % equivalent to UX*Mix1*UX
  end
  
elseif isequal(IncScheme,[1 -1 2])
  
  UtUX_ = conj(diagUX)*diagUX.';
  UY_ = diagUY*diagUY';
  % Pre-propagate mixing propagator to end of second period (= start of
  % experiment)
  MixX = diag(diagUX.^n(1))*Mix{1};
  Mix2 = Mix{2};
  for ix = 1:n(1)
    M = Mix2*MixX;
    FinalDensity = M*Density*M';
    for iy = 1:n(2)
      if nDetectors==1
        Signal(ix,iy) = Detector*FinalDensity(:);
      else
        for iDet = 1:nDetectors
          Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
        end
      end
      FinalDensity = UY_.*FinalDensity; % equivalent to UY*FinalDensity*UY'
    end
    MixX = UtUX_.*MixX; % equivalent to U^-1*MixX*U
  end
  
elseif isequal(IncScheme,[1 2 1])
  
  Mix1 = Mix{1};
  Mix2 = Mix{2};
  UUX_ = diagUX*diagUX.';
  UY = diag(diagUY);
  for iy = 1:n(2)
    MixY = Mix2*Mix1;
    MixYadj = MixY';
    for ix = 1:n(1)
      FinalDensity = MixY*Density*MixYadj;
      if nDetectors==1
        Signal(ix,iy) = Detector*FinalDensity(:);
      else
        for iDet = 1:nDetectors
          Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
        end
      end
      MixY = UUX_.*MixY; % equivalent to UX*MixY*UX
    end
    Mix1 = UY*Mix1;
  end
  
elseif isequal(IncScheme,[1 1 2 2])
  
  UUX_ = diagUX*diagUX.';
  UUY_ = diagUY*diagUY.';
  Mix2 = Mix{2};
  Uy3y = Mix{3};
  for iy = 1:n(2)
    Uy3y2 = Uy3y*Mix2;
    Ux1x = Mix{1};
    for ix = 1:n(1)
      P = Uy3y2*Ux1x;
      FinalDensity = P*Density*P';
      if nDetectors==1
        Signal(ix,iy) = Detector*FinalDensity(:);
      else
        for iDet = 1:nDetectors
          Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
        end
      end
      Ux1x = UUX_.*Ux1x; % equivalent to UX*Ux1x*UX
    end
    Uy3y = UUY_.*Uy3y; % equivalent to UY*Uy3y*UY
  end
  
elseif isequal(IncScheme,[1 2 2 1])
  
  Mix1 = Mix{1};
  Mix2 = Mix{2};
  Mix3 = Mix{3};
  UUX_ = diagUX*diagUX.';
  UUY_ = diagUY*diagUY.';
  for iy = 1:n(2)
    MixY = Mix3*Mix2*Mix1;
    MixYadj = MixY';
    for ix = 1:n(1)
      FinalDensity = MixY*Density*MixYadj;
      if nDetectors==1
        Signal(ix,iy) = Detector*FinalDensity(:);
      else
        for iDet = 1:nDetectors
          Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
        end
      end
      MixY = UUX_.*MixY; % equivalent to UX*MixY*UX
    end
    Mix2 = UUY_.*Mix2; % equivalent to UY*Mix2*UY
  end
  
elseif isequal(IncScheme,[1 2 -2 1])
  
  Mix1 = Mix{1};
  Mix2 = diag(diagUY.^n(2))*Mix{2}; % pre-propagate to endpoint of third delay
  Mix3 = Mix{3};
  UUX_ = diagUX*diagUX.';
  UtUY_ = conj(diagUY)*diagUY.';
  for iy = 1:n(2)
    MixY = Mix3*Mix2*Mix1;
    MixYadj = MixY';
    for ix = 1:n(1)
      FinalDensity = MixY*Density*MixYadj;
      if nDetectors==1
        Signal(ix,iy) = Detector*FinalDensity(:);
      else
        for iDet = 1:nDetectors
          Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
        end
      end
      MixY = UUX_.*MixY; % equivalent to UX*MixY*UX
    end
    Mix2 = UtUY_.*Mix2; % equivalent to UY'*Mix2*UY
  end
  
elseif isequal(IncScheme,[1 -1 1 -1])
  
  Mix1X = diag(diagU.^n)*Mix{1}; % pre-propagate to endpoint of second delay
  Mix2 = Mix{2};
  Mix3X = diag(diagU.^n)*Mix{3}; % pre-propagate to endpoint of fourth delay
  UtU_ = conj(diagU)*diagU.'; % propagator for Mix1 and Mix3 (add before, remove after)
  for ix = 1:n
    MixX = Mix3X*Mix2*Mix1X;
    FinalDensity = MixX*Density*MixX';
    if nDetectors==1
      Signal(ix) = Detector*FinalDensity(:);
    else
      for iDet = 1:nDetectors
        Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
      end
    end
    Mix1X = UtU_.*Mix1X; % equivalent to U'*Mix1X*U
    Mix3X = UtU_.*Mix3X; % equivalent to U'*Mix3X*U
  end
  
elseif isequal(IncScheme,[1 1 -1 -1])
  
  Mix1 = Mix{1};
  Mix2X = diag(diagU.^n)*Mix{2}; % pre-propagate to endpoint of third delay
  Mix3X = diag(diagU.^n)*Mix{3}; % pre-propagate to endpoint of fourth delay
  UU1_ = diagU*diagU.'; % propagator for Mix1 (add before and after)
  UU3_ = conj(diagU*diagU.'); % propagator for Mix3 (remove before and after)
  for ix = 1:n
    MixX = Mix3X*Mix2X*Mix1;
    FinalDensity = MixX*Density*MixX';
    if nDetectors==1
      Signal(ix) = Detector*FinalDensity(:);
    else
      for iDet = 1:nDetectors
        Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
      end
    end
    Mix1 = UU1_.*Mix1; % equivalent to U*Mix1*U
    Mix3X = UU3_.*Mix3X; % equivalent to U'*Mix3X*U'
  end
   
elseif isequal(IncScheme,[1 -1 -1 1])
  
  Mix1X = diag(diagU.^n)*Mix{1}; % pre-propagate to endpoint of second delay
  Mix2 = Mix{2};
  Mix3X = Mix{3}*diag(diagU.^n); % forward-propagate to endpoint of fourth delay
  UtU1_ = conj(diagU)*diagU.'; % propagator for Mix1 (add before, remove after)
  UtU3_ = diagU*diagU'; % propagator for Mix3 (remove before, add after)
  for ix = 1:n
    MixX = Mix3X*Mix2*Mix1X;
    FinalDensity = MixX*Density*MixX';
    if nDetectors==1
      Signal(ix) = Detector*FinalDensity(:);
    else
      for iDet = 1:nDetectors
        Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
      end
    end
    Mix1X = UtU1_.*Mix1X; % equivalent to U'*Mix1X*U
    Mix3X = UtU3_.*Mix3X; % equivalent to U*Mix3X*U'
  end
  
elseif isequal(IncScheme,[1 1 -1 -1 2])
  
  Mix1 = Mix{1};
  Mix2X = diag(diagUX.^n(1))*Mix{2}; % pre-propagate to endpoint of third delay
  Mix3X = diag(diagUX.^n(1))*Mix{3}; % pre-propagate to endpoint of fourth delay
  Mix4 = Mix{4};
  UU1_ = diagUX*diagUX.'; % propagator for Mix1 (add before and after)
  UU3_ = conj(diagUX*diagUX.'); % propagator for Mix3 (remove before and after)
  UY_ = diagUY*diagUY';
  for ix = 1:n(1)
    MixX = Mix4*Mix3X*Mix2X*Mix1;
    FinalDensity = MixX*Density*MixX';
    for iy = 1:n(2)
      if nDetectors==1
        Signal(ix,iy) = Detector*FinalDensity(:);
      else
        for iDet = 1:nDetectors
          Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
        end
      end
      FinalDensity = UY_.*FinalDensity; % equivalent to UY*FinalDensity*UY'
    end
    Mix1 = UU1_.*Mix1; % equivalent to U*Mix1*U
    Mix3X = UU3_.*Mix3X; % equivalent to U'*Mix3X*U'
  end
  
elseif isequal(IncScheme,[1 -1 -1 1 2])
  
  Mix1X = diag(diagUX.^n(1))*Mix{1}; % pre-propagate to endpoint of second delay
  Mix2 = Mix{2};
  Mix3X = Mix{3}*diag(diagUX.^n(1)); % forward-propagate to endpoint of fourth delay
  Mix4 = Mix{4};
  UtU1_ = conj(diagUX)*diagUX.'; % propagator for Mix1 (add before, remove after)
  UtU3_ = diagUX*diagUX'; % propagator for Mix3 (remove before, add after)
  UY_ = diagUY*diagUY';
  for ix = 1:n(1)
    MixX = Mix4*Mix3X*Mix2*Mix1X;
    FinalDensity = MixX*Density*MixX';
    for iy = 1:n(2)
      if nDetectors==1
        Signal(ix,iy) = Detector*FinalDensity(:);
      else
        for iDet = 1:nDetectors
          Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
        end
      end
      FinalDensity = UY_.*FinalDensity; % equivalent to UY*FinalDensity*UY'
    end
    Mix1X = UtU1_.*Mix1X; % equivalent to U'*Mix1X*U
    Mix3X = UtU3_.*Mix3X; % equivalent to U*Mix3X*U'
  end
  
else
  
  IncSchemeString = sprintf('%d ',IncScheme);
  error('The incrementation scheme\n  [%s]\ngiven in IncScheme is not supported.',IncSchemeString);
  
end

return
