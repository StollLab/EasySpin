function [Sys,Sigma0,DetOps,Events,Relaxation] = s_propagationsetup(Sys,Events,Opt)
% Spin system setup and input checking for saffron/spidyan

% Check spin system structure
% --------------------------------------------------------------------------------------
% Check for isotope mixtures
if ~isfield(Sys,'Nucs'), Sys.Nucs = ''; end

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
hasRelaxation = (isfield(Sys,'T1') && any(Sys.T1(:))) || (isfield(Sys,'T2') && any(Sys.T2(:)));
if isfield(Sys,'eqState') && ~hasRelaxation
  warning('An equilibrium state was provided, but no relaxation times. Equilibrium state ignored.');
end

if isfield(Sys,'T1T2')
  error('Sys.T1T2 is not supported. Provide relaxation times in Sys.T1 and Sys.T2.');
end

% Check if T1, T2 are provided or switch off/set to a very large value
if ~isfield(Sys,'T1') && ~isfield(Sys,'T2')
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

% Build initial state
if isempty(Sys.initState)
  % all electrons are -Sz, nuclei are not defined
  Sigma0 = -totSpinOps{3};
else
  % A user-provided initial state can depend on the Hamiltonian (and therefore
  % on the orientation and the simulation frame), so the caller sets it up
  % via initstate. Sigma0 is empty here, and so is the default equilibrium state.
  Sigma0 = [];
end

% Set up equilibrium state for relaxation (the relaxation superoperator
% depends on the Hamiltonian and is set up by the caller)

% --------------------------------------------------------------------------------------
if all([Sys.T1(:); Sys.T2(:)] == 0)
  Relaxation = [];
else
  logmsg(1,'  setting up the equilibrium state');

  if ~isscalar(Sys.T1) && ~isequal(size(Sys.T1),[Sys.nStates Sys.nStates])
    error('Sys.T1 must be a single number or a %dx%d matrix.',Sys.nStates,Sys.nStates);
  end
  if ~isscalar(Sys.T2) && ~isequal(size(Sys.T2),[Sys.nStates Sys.nStates])
    error('Sys.T2 must be a single number or a %dx%d matrix.',Sys.nStates,Sys.nStates);
  end

  Relaxation.Gamma = [];

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
    % equilibrium state is copied from initial state (empty if the caller sets
    % up the initial state, in which case the caller also sets this)
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
% The excitation operator X of each pulse is the total electron Sx, or a
% custom Hermitian operator (sop string or matrix). It couples to the
% x(Lab) component of the field. With circular polarization, the field
% also has a component along y(Lab), which couples to Y, the excitation
% operator rotated by 90 degrees about z(Lab): Y = Rz*X*Rz' with
% Rz = exp(-1i*pi/2*Fz), where Fz is the total Sz of all spins. For Sx,
% Y is Sy.
%
% The pulse Hamiltonian is H1 = Re(c)*xOp + Im(c)*yOp, where c is the
% pulse waveform (IQ) in the simulation frame, and yOp is empty if H1 has
% a single operator:
% - Lab frame: H1 = 2*Re(c)*X for linear polarization (field 2*Re(c)
%   along x(Lab)), and H1 = Re(c)*X + Im(c)*Y for circular polarization.
% - Shifted frame (rotating at f_s about the total electron Sz): the
%   elements of X with electron coherence order M = mS(i)-mS(j) acquire
%   exp(1i*M*2*pi*f_s*t). Only the part of the field that rotates with the
%   frame (rotating-wave approximation, as for the secular Hamiltonian) is
%   kept: H1 = conj(c)*Xp + c*Xp', where Xp is the M = +1 part of X
%   (linear polarization) or of (X+1i*Y)/2 (circular polarization). For
%   electron operators such as Sx, both are the same, and H1 =
%   Re(c)*Sx + Im(c)*Sy. Elements with M ~= +-1 oscillate at about the
%   microwave frequency or faster and are not driven by a microwave pulse.
%   They are rejected.
%   With c = a*exp(1i*psi), H1 = a*xOp.*exp(-1i*psi*M) is xOp rotated about
%   the total electron Sz, which commutes with the secular Hamiltonian.
%   M is stored in the event (CoherenceOrder), so that s_thyme can obtain
%   the propagators from a lookup table over the amplitude a.
% -------------------------------------------------------------------------

% Rz is diagonal, so the rotation multiplies each element by a phase
mTotal = 0;
for iSpin = 1:numel(Sys.Spins)
  mTotal = mTotal + full(diag(sop(Sys.Spins,[iSpin,3])));
end
rotPhases = exp(-1i*pi/2*(mTotal-mTotal.'));

% Electron coherence order of each matrix element
mS = full(diag(totSpinOps{3}));
M = mS-mS.';

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
        xOp = sop(Sys.Spins,Opt.ExcOperator{iPulse});
      elseif any(size(Opt.ExcOperator{iPulse}) ~= size(totSpinOps{1}))
        message = ['The excitation operator that you provided for pulse no. ' num2str(iPulse) ' does not have the same size as the density matrix.'];
        error(message)
      else
        xOp = Opt.ExcOperator{iPulse};
      end
      if norm(xOp-xOp','fro') > 1e-10*norm(xOp,'fro')
        error('The excitation operator for pulse no. %d is not Hermitian.',iPulse);
      end
      xOp = (xOp+xOp')/2;
    else
      % ... if not, the total Sx operator is used
      xOp = totSpinOps{1};
    end
    circular = strcmp(Events{iEvent}.Polarization,'circular');
    if circular
      yOp = xOp.*rotPhases;
    end
    Events{iEvent}.CoherenceOrder = [];
    if Events{iEvent}.FrameShift==0
      if circular
        Events{iEvent}.xOp = xOp;
        Events{iEvent}.yOp = yOp;
      else
        Events{iEvent}.xOp = 2*xOp;
        Events{iEvent}.yOp = [];
      end
    else
      if norm(xOp(abs(M)~=1)) > 1e-10*norm(xOp,'fro')
        error(['The excitation operator for pulse no. %d has elements that ' ...
          'do not change the electron mS by +-1. A microwave pulse cannot drive ' ...
          'them in a shifted simulation frame. Use Opt.SimFreq = 0 to ' ...
          'simulate in the lab frame.'],iPulse);
      end
      if circular
        Xp = (xOp+1i*yOp)/2;
      else
        Xp = xOp;
      end
      Xp(M~=1) = 0;
      Events{iEvent}.xOp = Xp+Xp';
      Events{iEvent}.yOp = 1i*(Xp'-Xp);
      Events{iEvent}.CoherenceOrder = M;
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