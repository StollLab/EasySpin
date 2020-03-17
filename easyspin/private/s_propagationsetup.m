function [Sys,Sigma0,DetOps,Events,Relaxation] = s_propagationsetup(Sys,Events,Opt)
% Spin system setup and input checking for saffron/spidyan
%

% Check spin system structure
% --------------------------------------------------------------------------------------
% Check for isotope mixtures
if ~isfield(Sys,'Nucs'), Sys.Nucs = ''; end
isoList = isotopologues(Sys.Nucs);
if numel(isoList)>1
  error('saffron/spidyan do not support isotope mixtures. Please specify pure isotopes in Sys.Nucs.');
end

% Validate spin system
[Sys,err] = validatespinsys(Sys);
error(err);
if Sys.MO_present, error('saffron/spidyan do not support general parameters!'); end
if any(Sys.L(:)), error('saffron/spidyan do not support L!'); end

if (Sys.nNuclei==0)
  logmsg(1,'  spin system: %d electron(s)',Sys.nElectrons);
else
  logmsg(1,'  spin system: %d electron(s), %d nuclei',Sys.nElectrons,Sys.nNuclei);
end

if isfield(Sys,'n') && any(Sys.n~=1)
  error('saffron/spidyan do not support Sys.n. Specify equivalent nuclei explicitly.')
end

if isfield(Sys,'nn') && any(Sys.nn(:)~=0)
  error('saffron/spidyan do not support nuclear-nuclear couplings (Sys.nn).');
end

logmsg(1,'  looking for and verifying fields related to relaxation');
% Look for equilibrium state and issue a warning if an equilibrium state is
% given, but no relaxation times
if isfield(Sys,'eqState') && ~isfield(Sys,'T2') && ~isfield(Sys,'T1')
  warning('An equilibrium state was provided, but no relaxation times. Equilibrium state ignored.');
end

% Check if T1, T2 are provided or switch off/set to a very large value
if ~isfield(Sys,'T1') && ~isfield(Sys,'T2')
  if isfield(Sys,'T1T2')
    error('T1 and T2 should be provided in Sys.T1 and Sys.T2 and not in Sys.T1T2.');
  end
  Sys.T1 = 0;
  Sys.T2 = 0;
elseif isfield(Sys,'T1') && ~isfield(Sys,'T2')
  Sys.T2 = 0;
elseif isfield(Sys,'T2') && ~isfield(Sys,'T1')
    Sys.T1 = 0;
end

if any([Sys.T1(:); Sys.T2(:)]<0) || any(~isreal([Sys.T1(:); Sys.T2(:)]))
  error('T1 and T2 in Sys.T1 and Sys.T2 must be positive and in microseconds.');
end

% Set up spin operators
% --------------------------------------------------------------------------------------
% Total electron spin operators, used for building initial and equilibrium
% denisity matrix as well es detection and excitation operators
logmsg(1,'  obtaining spin operators');
totSpinOps = cell(1,3);
for iSpin = 1:numel(Sys.S)
  if iSpin == 1
    for i = 1:3
      totSpinOps{i} = sop(Sys.Spins,[iSpin,i]);
    end
  else
    for i = 1:3
      totSpinOps{i} = totSpinOps{i} + sop(Sys.Spins,[iSpin,i]);
    end
  end
end

% Set up initial state
% --------------------------------------------------------------------------------------
logmsg(1,'  setting up the initial state');

% Build or load initial state
if isfield(Sys,'initState') && ~isempty(Sys.initState)
  % if some initial state was provided, this checks if the dimensions are
  % correct
  [a, b] = size(Sys.initState);
  if ischar(Sys.initState)
    error('String input for initial state not yet supported.')
  elseif Sys.nStates ~= a || Sys.nStates ~= b
    error('Initial state has to be a density matrix.')
  end
  Sigma0 = Sys.initState;
else
  % builds initial state, all electrons are -Sz, nuclei are not defined
  Sigma0 = -totSpinOps{3};
end

% Set up relaxation superoperator and equilibrium state

% --------------------------------------------------------------------------------------
if all([Sys.T1(:); Sys.T2(:)] == 0)
  Relaxation = [];
else
  logmsg(1,'  setting up the relaxation superoperator and the equilibrium state');

  Relaxation.Gamma = s_relaxationsuperoperator(Sys);
  
  % Set up equilibrium state (required for relaxation)
  if isfield(Sys,'eqState') && ~isempty(Sys.eqState)
    % if eqilibrium state was provided, this checks if the dimensions are correct
    [a, b] = size(Sys.eqState);
    if ischar(Sys.eqState)
      error('String input for equilibrium state not yet supported.')
    elseif Sys.nStates ~= a || Sys.nStates ~= b
      error('User-provided equilibrium state does not have the same dimensions as the density matrix.')
    end
    Relaxation.equilibriumState = Sys.eqState;
  else
    % initial state is copied from initial state
    Relaxation.equilibriumState  = Sigma0;
  end
end

if isfield(Opt,'Relaxation') && any(Opt.Relaxation) && isempty(Relaxation)
  error('You need to provide relaxation times Sys.T1 and Sys.T2 if you request relaxation.')
end

% Set up excitation operators
% --------------------------------------------------------------------------------------
logmsg(1,'  setting up the excitation operators');

% -------------------------------------------------------------------------
% Build the excitation operator for each pulse - if a custom excitation
% operator is provided in string form, sop will be called. If a matrix is
% provided, this matrix is used as excitation operator. If no custom
% excitation operators are given and no complex excitation was requested, 
% Sx is being used, for complex excitation Sx + Sy
% -------------------------------------------------------------------------

% initialize pulse counting
iPulse = 1;

% Loop over events and check if they are pulses
for iEvent = 1: length(Events)
  if strcmp(Events{iEvent}.type,'pulse')
    % Checks if user defined excitation operators were provided....
    if isfield(Opt,'ExcOperator') && ~isempty(Opt.ExcOperator) && (iPulse <= length(Opt.ExcOperator)) && ~isempty(Opt.ExcOperator{iPulse})
      % ... if yes, they are translated/verified and stored into the
      % current event ...
      if ischar(Opt.ExcOperator{iPulse})
        Events{iEvent}.xOp = sop(Sys.Spins,Opt.ExcOperator{iPulse});
      elseif any(size(Opt.ExcOperator{iPulse}) ~= size(Sigma0))
        message = ['The excitation operator that you provided for pulse no. ' num2str(iPulse) ' does not have the same size as the density matrix.'];
        error(message)
      else
         Events{iEvent}.xOp = Opt.ExcOperator{iPulse};
      end
    else
      % ... if not, the total Sx operator is used ...
      Events{iEvent}.xOp = totSpinOps{1};
      % ... and if complex excitation is required, the Sy operator is added
      % on top
      if Events{iEvent}.ComplexExcitation
        Events{iEvent}.xOp = Events{iEvent}.xOp + totSpinOps{2};
      end
    end
    % Increment pulse counter
    iPulse = iPulse + 1;
  end
end

% Set up detection operators
% -------------------------------------------------------------------------
logmsg(1,'  setting up the detection operators');

useDefaultDetOperator = ~isfield(Opt,'DetOperator') || isempty(Opt.DetOperator);


if useDefaultDetOperator
  % default value, equals to S+ for all electron spins 
  DetOps{1} = totSpinOps{1}+1i*totSpinOps{2}; 
else
  
  % catch the case Exp.DetOperator = 'z1' instead of {'z1'} - spidyan
  % checks and corrects for this case before
  if ischar(Opt.DetOperator)
    Opt.DetOperator = {Opt.DetOperator};
  end

  nDetOps = length(Opt.DetOperator);
  
  % compute DetOperators if input was string or validate them if input was
  % matrix
  DetOps = cell(1,nDetOps);  
  for iDetectionOperator = 1 : nDetOps
    if ischar(Opt.DetOperator{iDetectionOperator})
      DetOps{iDetectionOperator} = sop(Sys.Spins,Opt.DetOperator{iDetectionOperator});
    else
      if size(Opt.DetOperator{iDetectionOperator})~=size(totSpinOps{1})
        error('The size of the provided detection operator does not match the defined spin system.');
      end
      DetOps{iDetectionOperator} = Opt.DetOperator{iDetectionOperator};
    end
  end
logmsg(1,'  spin system validation finished successfully!');
end