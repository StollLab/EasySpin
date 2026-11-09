% sf_propagationsetup  Spin system setup and input checking for the full-space engine
%
%   [Sys,Sigma0,Det,Events,Relaxation] = sf_propagationsetup(Sys,Events,Opt)
%
%   Validates the spin system and sets up the initial state, the excitation
%   operators of the pulses, the detection operators and the equilibrium
%   state for sf_propagate (full-space engine of saffron).
%
%   Input:
%     Sys     spin system, including optional T1, T2 (µs), eqState, initState
%     Events  cell array of events, from sf_buildevents
%     Opt     options structure, with fields
%       DetOperator  detection operators (cell array of sop strings or
%                    matrices), default S+ of all electrons
%       DetFreq      detection frame frequency (GHz)
%       DetPhase     detection phase (radians) of each detection operator
%       Relaxation   relaxation on/off
%
%   Output:
%     Sys         validated spin system, with T1 and T2 (0 if not given)
%     Sigma0      initial density matrix, -Sz of all electrons; empty for
%                 Sys.initState, then set up by the caller
%     Det         detection structure, with fields
%       Ops          detection operators, cell array
%       M            electron coherence order of each matrix element
%       Phase        detection phase (radians) of each detection operator
%       FrameOffset  detection frame frequency relative to the simulation
%                    frame (GHz)
%     Events      input events, for pulses with the additional fields
%       xOp, yOp        excitation operators; yOp is empty if the pulse
%                       Hamiltonian has a single operator
%       CoherenceOrder  electron coherence orders of the elements of xOp in
%                       a shifted frame, otherwise empty
%       IdealOps        ideal pulses: total electron Sx and Sy, cell array
%     Relaxation  relaxation structure, [] without relaxation, with fields
%       Gamma             empty, set up by the caller
%       equilibriumState  equilibrium density matrix (empty if set up by
%                         the caller)

function [Sys,Sigma0,Det,Events,Relaxation] = sf_propagationsetup(Sys,Events,Opt)

% Check spin system structure
% --------------------------------------------------------------------------------------
% Check for isotope mixtures
if ~isfield(Sys,'Nucs'), Sys.Nucs = ''; end

% Validate spin system
[Sys,err] = validatespinsys(Sys);
error(err);
if Sys.MO_present, error('saffron does not support general parameters!'); end
if any(Sys.L(:)), error('saffron does not support L!'); end

if (Sys.nNuclei==0)
  logmsg(1,'  spin system: %d electron(s)',Sys.nElectrons);
else
  logmsg(1,'  spin system: %d electron(s), %d nuclei',Sys.nElectrons,Sys.nNuclei);
end

if isfield(Sys,'n') && any(Sys.n~=1)
  error('saffron does not support Sys.n. Specify equivalent nuclei explicitly.')
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
% The excitation operator X of each microwave pulse is the total electron
% Sx. It couples to the x(Lab) component of the field. With circular
% polarization, the field also has a component along y(Lab), which couples
% to Y, the excitation
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
%   Sx, both are the same, and H1 = Re(c)*Sx + Im(c)*Sy.
%   With c = a*exp(1i*psi), H1 = a*xOp.*exp(-1i*psi*M) is xOp rotated about
%   the total electron Sz, which commutes with the secular Hamiltonian.
%   M is stored in the event (CoherenceOrder), so that sf_propagate can obtain
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

% RF excitation operator (MHz/mT): nuclear Zeeman interaction with an RF
% field along x(Lab), H1 = 2*Re(IQ)*X0 with the lab-frame RF waveform IQ
% (mT). It has electron coherence order 0, so it is the same in shifted
% frames. The hyperfine enhancement via the electron Zeeman interaction
% with the RF field is not included.
X0 = 0;
for iNuc = 1:Sys.nNuclei
  X0 = X0 - Sys.gn(iNuc)*nmagn/1e3/planck/1e6*sop(Sys.Spins,[Sys.nElectrons+iNuc,1]);
end
X0 = full(X0);

% Excitation operators of the pulses
for iEvent = 1: length(Events)
  if strcmp(Events{iEvent}.type,'pulse') && isfield(Events{iEvent},'Ideal') && Events{iEvent}.Ideal
    % Ideal pulses rotate all electron spins (see sf_propagate)
    Events{iEvent}.IdealOps = {full(totSpinOps{1}), full(totSpinOps{2})};
    Events{iEvent}.xOp = [];
    Events{iEvent}.yOp = [];
    Events{iEvent}.CoherenceOrder = [];
    continue
  end
  if strcmp(Events{iEvent}.type,'pulse') && isfield(Events{iEvent},'Channel') && strcmp(Events{iEvent}.Channel,'rf')
    if Sys.nNuclei==0
      error('RF pulses need nuclei in the spin system.');
    end
    Events{iEvent}.xOp = 2*X0;
    Events{iEvent}.yOp = [];
    Events{iEvent}.CoherenceOrder = [];
    continue
  end
  if strcmp(Events{iEvent}.type,'pulse')
    xOp = totSpinOps{1};  % total electron Sx
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
  % Opt.DetOperator is a cell array (saffron wraps a single operator)
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
end

% Detection frame: signals and states are returned in the frame rotating
% at Opt.DetFreq (see sf_propagate). Det.Phase is the phase of the detection
% reference for each detection operator (Exp.DetPhase).
Det.Ops = DetOps;
Det.M = M;
if isfield(Opt,'DetPhase') && ~isempty(Opt.DetPhase)
  Det.Phase = Opt.DetPhase;
else
  Det.Phase = zeros(1,numel(DetOps));
end
if ~isscalar(Opt.DetFreq) || ~isreal(Opt.DetFreq)
  error(['Exp.DetFreq must be a single frequency (in GHz). All signals and ' ...
    'density matrices are given in one detection frame; operators without ' ...
    'electron coherence, such as Sz, are not affected by it.']);
end
if Opt.DetFreq<0
  error('Exp.DetFreq cannot be negative.');
end
Det.FrameOffset = Opt.DetFreq - Events{1}.FrameShift;  % GHz
logmsg(1,'  spin system validation finished successfully!');