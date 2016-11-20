% chili    Simulation of cw EPR spectra in the slow motional regime
%
%   chili(Sys,Exp,Opt)
%   spc = chili(...)
%   [B,spc] = chili(...)
%
%   Computes the slow-motion cw EPR spectrum of systems with
%   one electron and one nuclear spin.
%
%   Sys: spin system structure
%
%     Sys.tcorr           rotational correlation time (in seconds)
%     Sys.logtcorr        log10 of rotational correlation time (in seconds)
%     Sys.Diff            diffusion rate (s^-1)
%     Sys.logDiff         log10 of diffusion rate (s^-1)
%
%         All fields can have 1 (isotropic), 2 (axial) or 3 (rhombic) elements.
%         Precedence: logtcorr > tcorr > logDiff > Diff.
%
%     Sys.DiffFrame       Euler angles of the diffusion tensor (default [0 0 0])
%     Sys.lw              vector with FWHM residual broadenings
%                         1 element:  GaussianFWHM
%                         2 elements: [GaussianFWHM LorentzianFWHM]
%                         field sweep: mT, frequency sweep: MHz
%     Sys.lwpp            peak-to-peak line widths, same format as Sys.lw
%     Sys.lambda          ordering potential coefficients
%                         [lambda20 lambda22 lambda40 lambda42 lambda44]
%     Sys.Exchange        Heisenberg exchange frequency (MHz)
%
%    Exp: experimental parameter settings
%      mwFreq              microwave frequency, in GHz (for field sweeps)
%      Range               sweep range, [sweepmin sweepmax], in mT (for field sweep)
%      CenterSweep         sweep range, [center sweep], in mT (for field sweeps
%      Field               static field, in mT (for frequency sweeps)
%      mwRange             sweep range, [sweepmin sweepmax], in GHz (for freq. sweeps)
%      mwCenterSweep       sweep range, [center sweep], in GHz (for freq. sweeps)
%      nPoints             number of points
%      Harmonic            detection harmonic: 0, 1 (default), 2
%      ModAmp              peak-to-peak modulation amplitude, in mT (field sweeps only)
%      mwPhase             detection phase (0 = absorption, pi/2 = dispersion)
%      Temperature         temperature, in K
%
%   Opt: simulation options
%      LLKM            basis size: [evenLmax oddLmax Kmax Mmax]
%      PostConvNucs    nuclei to include perturbationally via post-convolution
%      Verbosity       0: no display, 1: show info
%      nKnots          number of knots for powder simulation
%      Symmetry        symmetry to use for powder simulation
%
%   Output:
%     B      magnetic field axis vector, in mT
%     spc    simulated spectrum, arbitrary units
%
%     If no output arguments are specified, chili plots the
%     simulated spectrum.

function varargout = chili(Sys,Exp,Opt)


if (nargin==0), help(mfilename); return; end

error(chkmlver);
if (nargin<2) || (nargin>3), error('Wrong number of input arguments!'); end
if (nargout<0), error('Not enough output arguments.'); end
if (nargout>2), error('Too many output arguments.'); end

if (nargin<3), Opt = struct('unused',NaN); end

if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 0; % Log level
end

% --------License ------------------------------------------------
LicErr = 'Could not determine license.';
Link = 'epr@eth'; eschecker; error(LicErr); clear Link LicErr
% --------License ------------------------------------------------

global EasySpinLogLevel;
EasySpinLogLevel = Opt.Verbosity;

%==================================================================
% Loop over components and isotopologues
%==================================================================
FrequencySweep = ~isfield(Exp,'mwFreq') & isfield(Exp,'Field');

if FrequencySweep
  SweepAutoRange = (~isfield(Exp,'mwRange') || isempty(Exp.mwRange)) && ...
    (~isfield(Exp,'mwCenterSweep') || isempty(Exp.mwCenterSweep));
else
  SweepAutoRange = (~isfield(Exp,'Range') || isempty(Exp.Range)) && ...
    (~isfield(Exp,'CenterSweep') || isempty(Exp.CenterSweep));
end
if ~isfield(Opt,'IsoCutoff'), Opt.IsoCutoff = 1e-4; end

if ~isfield(Sys,'singleiso')
  
  [SysList,weight] = expandcomponents(Sys,Opt.IsoCutoff);
  
  if (numel(SysList)>1) && SweepAutoRange
    if FrequencySweep
      error('Multiple components: Please specify frequency range manually using Exp.mwRange or Exp.mwCenterSweep.');
    else
      error('Multiple components: Please specify field range manually using Exp.Range or Exp.CenterSweep.');
    end
  end
  
  spec = 0;
  for iComponent = 1:numel(SysList)
    [xAxis,spec_] = chili(SysList{iComponent},Exp,Opt);
    spec = spec + spec_*weight(iComponent);
  end
  
  % Output and plotting
  switch nargout
    case 0
      cla
      if FrequencySweep
        if (xAxis(end)<1)
          plot(xAxis*1e3,spec);
          xlabel('frequency (MHz)');
        else
          plot(xAxis,spec);
          xlabel('frequency (GHz)');
        end
        title(sprintf('%0.8g mT',Exp.Field));
      else
        if (xAxis(end)<10000)
          plot(xAxis,spec);
          xlabel('magnetic field (mT)');
        else
          plot(xAxis/1e3,spec);
          xlabel('magnetic field (T)');
        end
        title(sprintf('%0.8g GHz',Exp.mwFreq));
      end
      axis tight
      ylabel('intensity (arb.u.)');
    case 1, varargout = {spec};
    case 2, varargout = {xAxis,spec};
  end
  return
end
%==================================================================


logmsg(1,'-- slow motion regime simulation ----------------------------------');

% Spin system
%-------------------------------------------------------------------
out = isotopologues(Sys.Nucs);
if (out.nIso>1)
  error('chili does not support isotope mixtures. Please specify pure isotopes in Sys.Nucs.');
end

[Sys,err] = validatespinsys(Sys);
error(err);
if Sys.MO_present, error('salt does not support general parameters!'); end
if any(Sys.L(:)), error('salt does not support L!'); end

if Sys.fullg
  idx = 1:3;
  for iElectron = 1:Sys.nElectrons
    mean_g(iElectron) = mean(eig(Sys.g(idx,:)));
    idx = idx + 3;
  end
  mT2MHz = mt2mhz(1,mean(mean_g));
else
  mT2MHz = mt2mhz(1,mean(mean(Sys.g)));
end

if any(Sys.HStrain(:)) || any(Sys.gStrain(:)) || any(Sys.AStrain(:)) || any(Sys.DStrain(:))
  error('chili does not support strains (HStrain, gStrain, AStrain, DStrain). Please remove from spin system.');
end

% Convolution with Gaussian only. Lorentzian broadening is 
% included in the slow-motion simulation via T2.
ConvolutionBroadening = any(Sys.lw(1)>0);

% Dynamics and ordering potential
%-------------------------------------------------------------------
if isfield(Sys,'psi')
  error('Sys.psi is obsolete. Remove it from your code. See the documentation for details.');
end

if ~isfield(Sys,'DiffFrame'), Sys.DiffFrame = [0 0 0]; end
if ~isfield(Sys,'Exchange'), Sys.Exchange = 0; end
if ~isfield(Sys,'lambda'), Sys.lambda = []; end

if isfield(Sys,'tcorr'), Dynamics.tcorr = Sys.tcorr; end
if isfield(Sys,'Diff'), Dynamics.Diff = Sys.Diff; end
if isfield(Sys,'logtcorr'), Dynamics.logtcorr = Sys.logtcorr; end
if isfield(Sys,'logDiff'), Dynamics.logDiff = Sys.logDiff; end
if isfield(Sys,'lwpp'), Dynamics.lwpp = Sys.lwpp; end
if isfield(Sys,'lw'), Dynamics.lw = Sys.lw; end

Dynamics.Exchange = Sys.Exchange;
Potential.lambda = Sys.lambda;
usePotential = ~isempty(Potential.lambda) && ~all(Potential.lambda==0);

% Experimental settings
%-------------------------------------------------------------------
if isfield(Exp,'MOMD')
  error('Exp.MOMD is obsolete. Remove it from your code. See the documentation for details.');
end

if ~isfield(Exp,'nPoints'), Exp.nPoints = 1024; end
if ~isfield(Exp,'Harmonic'), Exp.Harmonic = []; end
if ~isfield(Exp,'mwPhase'), Exp.mwPhase = 0; end
if ~isfield(Exp,'Temperature'), Exp.Temperature = NaN; end
if ~isfield(Exp,'ModAmp'), Exp.ModAmp = 0; end
if ~isfield(Exp,'Mode'), Exp.Mode = 'perpendicular'; end
if ~isfield(Exp,'Ordering'), Exp.Ordering = []; end
if ~isfield(Exp,'CrystalOrientation'), Exp.CrystalOrientation = []; end

% Number of points
if any(~isreal(Exp.nPoints)) || numel(Exp.nPoints)>1 || (Exp.nPoints<2)
  error('Problem with Exp.nPoints. Needs to be a number >= 2.')
end

% Temperature
if ~isnan(Exp.Temperature)
  if (numel(Exp.Temperature)~=1) || isinf(Exp.Temperature) || (Exp.Temperature<0)
    error('Problem with Exp.Temperature. If given, Exp.Temperature must be a positive value.')
  end
end

% Microwave frequency
if ~isfield(Exp,'mwFreq')
  if ~isfield(Exp,'Field')
    error('Please supply either the microwave frequency in Exp.mwFreq (for field sweeps) or the magnetic field in Exp.Field (for frequency sweeps).');
  end
  FieldSweep = false;
else
  if isfield(Exp,'Field')
    error('Give either Exp.mwFreq (for a field sweep) or Exp.Field (for a frequency sweep), but not both.');
  end
  FieldSweep = true;
end
if FieldSweep
  if (numel(Exp.mwFreq)~=1) || any(Exp.mwFreq<=0) || ~isreal(Exp.mwFreq)
    error('Uninterpretable microwave frequency in Exp.mwFreq.');
  end
  logmsg(1,'  field sweep, mw frequency %0.8g GHz',Exp.mwFreq);
else
  if (numel(Exp.Field)~=1) || any(Exp.Field<=0) || ~isreal(Exp.Field)
    error('Uninterpretable magnetic field in Exp.Field.');
  end
  logmsg(1,'  frequency sweep, magnetic field %0.8g mT',Exp.Field);
end

% Sweep range (magnetic field, or frequency)
if FieldSweep
  if isfield(Exp,'CenterSweep')
    if isfield(Exp,'Range')
      logmsg(0,'Using Exp.CenterSweep and ignoring Exp.Range.');
    end
  else
    if isfield(Exp,'Range')
      Exp.CenterSweep = [mean(Exp.Range) diff(Exp.Range)];
    else
      if (Sys.nElectrons==1) && (Sys.S==1/2)
        logmsg(1,'  automatic determination of sweep range');
        Stretch = 1.25;
        I = nucspin(Sys.Nucs).';
        if numel(I)>0
          Amax = max(abs(Sys.A),[],2);
          hf = sum(I.*Amax)*1e6; % MHz -> Hz
        else
          hf = 0;
        end
        gmax = max(Sys.g(:));
        gmin = min(Sys.g(:));
        if FieldSweep
          minB = planck*(Exp.mwFreq*1e9 - hf)/bmagn/gmax/1e-3;
          maxB = planck*(Exp.mwFreq*1e9 + hf)/bmagn/gmin/1e-3;
          Exp.CenterSweep = [(maxB+minB)/2, Stretch*max(maxB-minB,5)];
        else
          minE = bmagn*Exp.Field*1e-3*gmin/planck - hf; % Hz
          maxE = bmagn*Exp.Field*1e-3*gmax/planck + hf; % Hz
          Exp.CenterSweep = [(maxE+minE)/2, Stretch*max(maxE-minE,10e6)]/1e9; % GHz
        end
      else
        error('Cannot automatically determine sweep range for this spin system.');
      end
    end
  end
else
  if isfield(Exp,'mwCenterSweep')
    if isfield(Exp,'mwRange')
      logmsg(0,'Using Exp.mwCenterSweep and ignoring Exp.mwRange.');
    end
  else
    if isfield(Exp,'mwRange')
      Exp.mwCenterSweep = [mean(Exp.mwRange) diff(Exp.mwRange)];
    else
      error('Either Exp.mwRange or Exp.mwCenterSweep need to be given.');
    end
  end
end

if FieldSweep
  CenterField = Exp.CenterSweep(1);
  Sweep = Exp.CenterSweep(2);
  Exp.Range = Exp.CenterSweep(1) + [-1 1]/2*Sweep;
  if any(Exp.Range<0) || diff(Exp.Range)<=0
    error('Invalid sweep range! Check Exp.CenterSweep or Exp.Range.');
  end
else
  CenterFreq = Exp.mwCenterSweep(1);
  Sweep = Exp.mwCenterSweep(2);
  Exp.mwRange = Exp.mwCenterSweep(1) + [-1 1]/2*Sweep;
  CenterField = Exp.Field;
  if any(Exp.mwRange<0) || diff(Exp.mwRange)<=0
    error('Invalid sweep range! Check Exp.mwCenterSweep or Exp.mwRange.');
  end
end

if FieldSweep
  logmsg(1,'  field range (mT): min %g, max %g, center %g, width %g',...
    Exp.Range(1),Exp.Range(2),CenterField,Sweep);
else
  logmsg(1,'  frequency range (GHz): min %g, max %g, center %g, width %g',...
    Exp.mwRange(1),Exp.mwRange(2),CenterFreq,Sweep);
end

% Detection harmonic
if ~isfield(Exp,'Harmonic') || isempty(Exp.Harmonic) || isnan(Exp.Harmonic)
  if FieldSweep
    Exp.Harmonic = 1;
  else
    Exp.Harmonic = 0;
  end
end
if ~any(Exp.Harmonic==[-1,0,1,2])
  error('Exp.Harmonic must be 0, 1 or 2.');
end

% Modulation amplitude
if any(Exp.ModAmp<0) || any(isnan(Exp.ModAmp)) || numel(Exp.ModAmp)~=1
  error('Exp.ModAmp must be either a single positive number or zero.');
end
if (Exp.ModAmp>0)
  if FieldSweep
    logmsg(1,'  field modulation, amplitude %g mT',Exp.ModAmp);
    if (Exp.Harmonic<1)
      error('With field modulation (Exp.ModAmp), Exp.Harmonic=0 does not work.');
    end
    Exp.ModHarmonic = Exp.Harmonic;
    Exp.ConvHarmonic = 0;
    Exp.DerivHarmonic = 0;
  else
    error('Exp.ModAmp cannot be used with frequency sweeps.');
  end
else
  Exp.ModHarmonic = 0;
  if ConvolutionBroadening
    Exp.ConvHarmonic = Exp.Harmonic;
    Exp.DerivHarmonic = 0;
  else
    Exp.ConvHarmonic = 0;
    Exp.DerivHarmonic = Exp.Harmonic;
  end
end

% Resonator mode
switch Exp.Mode
  case 'perpendicular', ParallelMode = false;
  case 'parallel', ParallelMode = true;
  otherwise, error('Exp.Mode must be either ''perpendicular'' or ''parallel''.');
end
logmsg(1,'  harmonic %d, %s mode',Exp.Harmonic,Exp.Mode);

% Complain if fields only valid in pepper() are given
if isfield(Exp,'Orientations')
  warning('Exp.Orientations is obsolete. Use Exp.CrystalOrientations instead.');
end
if isfield(Exp,'CrystalSymmetry')
  warning('Exp.CrystalSymmetry is not used by chili.');
end

% Partial ordering
if ~isempty(Exp.Ordering)
  %if ~PowderSimulation
  %  error('Partial ordering (Exp.Ordering) can only be used in a powder simulation.');
  %end
  if isnumeric(Exp.Ordering) && (numel(Exp.Ordering)==1) && isreal(Exp.Ordering)
    UserSuppliedOrderingFcn = false;
    logmsg(1,'  partial order (built-in function, coefficient = %g)',Exp.Ordering);
  elseif isa(Exp.Ordering,'function_handle')
    UserSuppliedOrderingFcn = true;
    logmsg(1,'  partial order (user-supplied function)');
  else
    error('Exp.Ordering must be a single number or a function handle.');
  end
end

% Determine whether to do a powder simulation
if ~usePotential
  if isempty(Exp.Ordering) || all(Exp.Ordering==0)
    logmsg(1,'  No ordering potential given, skipping powder simulation.');
    PowderSimulation = false;
  else
  logmsg(1,'  Ordering potential given, doing powder simulation.');
    PowderSimulation = true;
  end    
else
  if ~isempty(Exp.CrystalOrientation)
    logmsg(1,'  Ordering potential given, doing single-crystal simulation.');
    PowderSimulation = false;
  else
    logmsg(1,'  Ordering potential given, doing powder simulation.');
    PowderSimulation = true;
  end
end

% Options
%-------------------------------------------------------------------
if isempty(Opt), Opt = struct('unused',NaN); end
if ~isfield(Opt,'Rescale'), Opt.Rescale = 1; end % rescale A before Lanczos
if ~isfield(Opt,'Threshold'), Opt.Threshold = 1e-6; end
if ~isfield(Opt,'Diagnostic'), Opt.Diagnostic = 0; end
if ~isfield(Opt,'Solver'), Opt.Solver = 'L'; end
if ~isfield(Opt,'Lentz'), Opt.Lentz = 1; end
if ~isfield(Opt,'IncludeNZI'), Opt.IncludeNZI = true; end
if ~isfield(Opt,'pqOrder'), Opt.pqOrder = false; end
if ~isfield(Opt,'Symmetry'), Opt.Symmetry = 'Dinfh'; end
if ~isfield(Opt,'SymmFrame'), Opt.SymmFrame = []; end
if ~isfield(Opt,'PostConvNucs'), Opt.PostConvNucs = ''; end

% Obsolete options
% Opt.MOMD was used prior to 5.0 for powder simulations (in the presence of ordering potential)
if isfield(Opt,'MOMD')
  error('Opt.MOMD is obsolete. Remove it from your code. Now, a powder simulation is automatically performed whenever an ordering potential is given.');
end

% Set default method for constructing Liouvillian
if ~isfield(Opt,'LiouvMethod') || isempty(Opt.LiouvMethod)
  if (Sys.nElectrons==1) && (Sys.S==1/2) && (Sys.nNuclei<=2)
    Opt.LiouvMethod = 'Freed';
  else
    Opt.LiouvMethod = 'general';
  end
else
end

[LiouvMethod,err] = parseoption(Opt,'LiouvMethod',{'Freed','general'});
error(err);
generalLiouvillian = (LiouvMethod==2);

if ~generalLiouvillian
  if (Sys.nElectrons>1) || (Sys.S~=1/2)
    error('Opt.LiouvMethod=''Freed'' does not work with this spin system.');
  end
end

if generalLiouvillian
  if usePotential
    error('Ordering potential not supported for Opt.LiouvMethod=''general''.');
  end
end

% Post-convolution nuclei
doPostConvolution = ~isempty(Opt.PostConvNucs);
if doPostConvolution
  Opt.PostConvNucs = sort(unique(Opt.PostConvNucs));
  if any(Opt.PostConvNucs<1) || any(Opt.PostConvNucs>Sys.nNuclei)
    error('Opt.PostConvNucs must contain indices of nuclei (1 to %d).',Sys.nNuclei);
  end
  nPostConvNucs = numel(Opt.PostConvNucs);
  if (Sys.nNuclei-nPostConvNucs>2) && ~generalLiouvillian
    error('Cannot have more than two nuclei for the Stochastic Liouville equation with this Opt.Method.');
  end
  fullSys = Sys;
  Sys = nucspinrmv(Sys,Opt.PostConvNucs);
  Sys.processed = 0;
  [Sys,err] = validatespinsys(Sys);
  error(err);
end

if ~generalLiouvillian
  % Pick functions for the calculation of the Liouvillian
  switch Sys.nNuclei
    case 0, chili_lm = @chili_lm0;
    case 1, chili_lm = @chili_lm1;
    case 2, chili_lm = @chili_lm2;
    otherwise
      error('Cannot have more than two nuclei for the Stochastic Liouville equation with this Opt.Method.');
  end
end

if any(Sys.n~=1)
  error('Cannot solve the Stochastic Liouville equation for systems with any Sys.n > 1.');
end

if ~isfield(Opt,'nKnots'), Opt.nKnots = [5 0]; end
if numel(Opt.nKnots)<1, Opt.nKnots(1) = 5; end
if numel(Opt.nKnots)<2, Opt.nKnots(2) = 0; end

% Basis settings
if isfield(Opt,'LLMK')
  error('Opt.LLMK is not a valid field. Use Opt.LLKM.');
end
if ~isfield(Opt,'LLKM')
  Opt.LLKM = [14 7 6 2];
end
Basis.LLKM = Opt.LLKM;
if ~isfield(Opt,'jKmin')
  Opt.jKmin = [];
end
Basis.jKmin = Opt.jKmin;
if ~isfield(Opt,'deltaK')
  Opt.deltaK = [];
end
Basis.deltaK = Opt.deltaK;

if ~isfield(Opt,'pSmin')
  Opt.pSmin = 0;
end
Basis.pSmin = Opt.pSmin;

% Maximum nuclear coherence order
if ~isfield(Opt,'pImax')
  Opt.pImax = [];
end
if Opt.pImax<0
  error('Opt.pImax must be 0 or larger.');
end
Basis.pImax = Opt.pImax;

if ~isfield(Opt,'MpSymm')
  Opt.MpSymm = false;
end
Basis.MpSymm = Opt.MpSymm;

switch Opt.Solver
  case 'L'
    if (Opt.Lentz==1) % Lentz method
      SolverString = 'Lanczos tridiagonalization, left-to-right continued fraction evaluation';
    else
      SolverString = 'Lanczos tridiagonalization, right-to-left continued fraction evaluation';
    end
  case 'C'
    SolverString = 'conjugate gradients tridiagonalization, right-to-left continued fraction evaluation';
  case 'R'
    SolverString = 'biconjugate gradients, stabilized';
  case '\'
    SolverString = 'backslash linear';
  case 'D'
    SolverString = 'direct method (eigenbasis, Binsch)';
  otherwise
    error('Unknown method in Options.Solver. Must be ''L'', ''R'', ''C'', or ''\''.');
end
logmsg(1,'Solver: %s',SolverString);

if ~generalLiouvillian
  maxElements = 5e6; % used in chili_lm
  maxRows = 2e5; % used in chili_lm
  if ~isfield(Opt,'Allocation')
    Opt.Allocation = [maxElements maxRows];
  elseif numel(Opt.Allocation)<2
    Opt.Allocation(2) = maxRows;
  end
  if Opt.Allocation(1)<1e3
    error('Opt.Allocation(1) (maximum number elements) is too small.');
  end
  logmsg(2,'  allocation: %d max elements, %d max rows',Opt.Allocation(1),Opt.Allocation(2));
end

% Process
%-------------------------------------------------------
if generalLiouvillian
  logmsg(1,'  using general Liouvillian code');
else
  logmsg(1,'  using S=1/2 Liouvillian code');
end
Sys = istospinsys(Sys,CenterField,Opt.IncludeNZI);

[Dynamics,err] = processdynamics(Dynamics,FieldSweep);
error(err);

% Ordering potential
%------------------------------------------------------------------
if ~isfield(Potential,'lambda'), Potential.lambda = [0 0 0 0 0]; end
if numel(Potential.lambda)<5, Potential.lambda(5) = 0; end
if numel(Potential.lambda)>5, error('Too many potential coefficients!'); end

Potential.L = [2 2 4 4 4];
Potential.K = [0 2 0 2 4];

% Calculate list of potential coefficients
Potential.xlk = chili_xlk(Potential,Dynamics.Diff);

% Basis
%------------------------------------------------------------------
Basis = processbasis(Sys,Basis,max(Potential.K));
if isempty(Basis.jKmin)
  error('Basis.jKmin is empty. Please report.');
end

% Set up horizontal sweep axis
% (nu is used internally, xAxis is used for user output)
if FieldSweep
  FreqSweep = Sweep*mT2MHz*1e6; % mT -> Hz
  nu = Exp.mwFreq*1e9 - linspace(-1,1,Exp.nPoints)*FreqSweep/2;  % Hz
  xAxis = linspace(Exp.Range(1),Exp.Range(2),Exp.nPoints);  % field axis, mT
else
  nu = linspace(Exp.mwRange(1),Exp.mwRange(2),Exp.nPoints)*1e9;  % Hz
  xAxis = nu/1e9; % frequency axis, GHz
end
omega0 = complex(1/(Dynamics.T2),2*pi*nu); % angular frequency

% Set up list of orientations
%=====================================================================
if (PowderSimulation)
  if Opt.nKnots(1)==1
    phi = 0;
    theta = 0;
    GridWeights = 4*pi;
  else
    [Vecs,GridWeights] = sphgrid(Opt.Symmetry,Opt.nKnots(1),'cf');
    % Transform vector to reference frame representation and convert to polar angles.
    if isempty(Opt.SymmFrame)
      [phi,theta] = vec2ang(Vecs);
    else
      [phi,theta] = vec2ang(Opt.SymmFrame*Vecs);
    end
  end
  logmsg(1,'  powder simulation with %d orientations',numel(phi));
else
  if ~isempty(Exp.CrystalOrientation)
    phi = Exp.CrystalOrientation(1);
    theta = Exp.CrystalOrientation(2);
  else
    phi = 0;
    theta = 0;
  end
  GridWeights = 4*pi;
  logmsg(2,'  single-orientation simulation');
end
nOrientations = numel(phi);
Basis.DirTilt = any(theta~=0);

% Partial ordering for protein/macromolecule
if ~isempty(Exp.Ordering)
  if (UserSuppliedOrderingFcn)
    OrderingWeights = feval(Exp.Ordering,phi,theta);
    if any(OrderingWeights)<0, error('User-supplied orientation distribution gives negative values!'); end
    if max(OrderingWeights)==0, error('User-supplied orientation distribution is all-zero.'); end
    logmsg(2,'  user-supplied ordering potential');
  else
    logmsg(2,'  standard ordering potential');
    U = -Exp.Ordering*plegendre(2,0,cos(theta)); % ordering potential
    OrderingWeights = exp(-U);
  end
else
  OrderingWeights = ones(1,nOrientations);
end

Weights = GridWeights.*OrderingWeights;
Weights = 4*pi*Weights/sum(Weights);


% Basis set preparations
%-----------------------------------------------------------------------
logmsg(1,'Setting up basis set...');
logmsg(1,'  spatial basis: Leven max %d, Lodd max %d, Kmax %d, Mmax %d, deltaK %d, jKmin %+d',...
  Basis.LLKM(1),Basis.LLKM(2),Basis.LLKM(3),Basis.LLKM(4),Basis.deltaK,Basis.jKmin);
logmsg(1,'  spin basis: pSmin %+d, pImax %d',Basis.pSmin,Basis.pImax);
logmsg(1,'  M-p symmetry: %d',Basis.MpSymm);

if generalLiouvillian
  
  % Set up basis
  Basis.List = generatebasis(Basis);
  nOriBasis = size(Basis.List,1);
  nSpinBasis = Sys.nStates^2;
  logmsg(1,'  complete product basis size: %d (%d spatial, %d spin)',nOriBasis*nSpinBasis,nOriBasis,nSpinBasis);
  
  % Index vector for reordering basis states from m1-m2 order (standard) to p-q order (Freed)
  [idxpq,mm,pq] = pqorder(Sys.Spins);
  
  % Removing unwanted spin functions
  rmv = false;
  % (1) remove any transitions with pS<pSmin
  for ie = 1:Sys.nElectrons
    rmv = rmv | pq(:,2*ie-1)<Basis.pSmin;
  end
  % (2) remove any transitions with |pI|>pImax, for each nucleus
  for in = 1:Sys.nNuclei
    rmv = rmv | any(abs(pq(:,2*Sys.nElectrons+2*in-1))>Basis.pImax(in),2);
  end
  keep = repmat(~rmv,nOriBasis,1);
  logmsg(1,'  pruning spin basis: keeping %d of %d functions',sum(~rmv),nSpinBasis);
  
  % Apply M=p-1 symmetry (Meirovitch Eq. (A47))
  if Opt.MpSymm
    M = Basis.List(:,2);
    psum = sum(pq(:,1:2:end),2);
    keep_Mp = bsxfun(@minus,psum,M.')==1; % keep only basis states with pS+pI-M == 1
    keep = keep & keep_Mp(:);
    logmsg(1,'  applying M-p symmetry: keeping %d of %d functions',sum(keep),numel(keep));
  end
  
  logmsg(1,'  final basis size: %d (%f%% of %d)',sum(keep),100*sum(keep)/nOriBasis/nSpinBasis,nOriBasis*nSpinBasis);
    
else
  
  [Basis.Size,Basis.SpatialSize,Indices] = chili_basiscount(Basis,Sys);
  logmsg(1,'  basis size: %d',Basis.Size);
    
end

% Precalculating operator matrices
%-----------------------------------------------------------------------
if generalLiouvillian
  
  logmsg(1,'Generating all cartesian spin operators');
  for iSpin = 1:numel(Sys.Spins)
    SpinOps{iSpin,1} = sop(Sys.Spins,iSpin,1,'sparse');
    SpinOps{iSpin,2} = sop(Sys.Spins,iSpin,2,'sparse');
    SpinOps{iSpin,3} = sop(Sys.Spins,iSpin,3,'sparse');
  end
  
  logmsg(1,'Generating ISTOs and precalculating 3j symbols');
  [T0,T1,T2,F0,F1,F2] = magint(Sys,SpinOps,CenterField,Opt.IncludeNZI);
  [jjj0,jjj1,jjj2] = jjjsymbol(Basis.LLKM,any(F1(:)));
  
  logmsg(1,'Setting up the detection operator');
  SxOp = SpinOps{1,1};
  for e = 2:Sys.nElectrons
    SxOp = SxOp + SpinOps{e,1};
  end
  SxOp = SxOp(:);
  if Opt.pqOrder
    SxOp = SxOp(idxpq);
  end
  
  logmsg(1,'Calculating the relaxation superoperator matrix');
  % Calculate relaxation superoperator in spatial basis, expand to full product
  % basis, and remove unwanted basis functions.
  Gamma = diffsuperop(Dynamics.Diff,Basis.List);
  Gamma = spkroneye(Gamma,Sys.nStates^2);
  Gamma = Gamma(keep,keep);
  
end

% Loop over all orientations
%=====================================================================
spec = 0;
for iOri = 1:nOrientations
  
  % Set up orientation
  %-------------------------------------------------------
  logmsg(2,'orientation %d of %d: phi = %g°, theta = %g° (weight %g)',...
    iOri,nOrientations,phi(iOri)*180/pi,theta(iOri)*180/pi,Weights(iOri));

  if generalLiouvillian
    D1 = wignerd(1,[phi(iOri),theta(iOri),0]);
    D2 = wignerd(2,[phi(iOri) theta(iOri) 0]);
    [Q0,Q1,Q2] = rbos(D1,D2,T0,T1,T2,F0,F1,F2);
    if Opt.pqOrder
      Q0 = Q0(idxpq,idxpq);
      for k = 1:numel(Q1)
        Q1{k} = Q1{k}(idxpq,idxpq);
      end
      for k = 1:numel(Q2)
        Q2{k} = Q2{k}(idxpq,idxpq);
      end
    end
  else
    Sys.d2psi = wignerd(2,[phi(iOri) theta(iOri) 0]);
  end
  
  % Starting vector
  %-------------------------------------------------------
  logmsg(1,'Computing starting vector...');
  if generalLiouvillian
    % set up in full product basis, then prune
    StartingVector = startvec(Basis.List,SxOp);
    StartingVector = StartingVector(keep);
  else
    StartingVector = chili_startingvector(Basis,Potential,Sys.I);
  end
  BasisSize = size(StartingVector,1);
  logmsg(1,'  vector size: %dx1',BasisSize);
  logmsg(1,'  non-zero elements: %d/%d (%0.2f%%)',...
    nnz(StartingVector),BasisSize,100*nnz(StartingVector)/BasisSize);
  logmsg(1,'  maxabs %g, norm %g',full(max(abs(StartingVector))),norm(StartingVector));
  
  % Liouvillian matrix
  %-------------------------------------------------------
  logmsg(1,'Computing Liouvillian matrix...');
  
  if generalLiouvillian
    H = liouvhamiltonian(Basis.List,Q0,Q1,Q2,jjj0,jjj1,jjj2);
    L = -2i*pi*H(keep,keep) + Gamma;
    nDim = size(L,1);
  else
    Sys.DirTilt = Basis.DirTilt; % used in chili_lm
    Dynamics.xlk = Potential.xlk; % used in chili_lm
    Dynamics.maxL = size(Potential.xlk,1)-1; % used in chili_lm
    [r,c,Vals,nDim,nElm] = chili_lm(Sys,Basis.v,Dynamics,Opt.Allocation);
    idx = 1:nElm;
    
    % Downsize oversized arrays to avoid running out of memory
    r = r(idx);
    c = c(idx);
    Vals = Vals(idx);
    
    L = sparse(r+1,c+1,Vals,BasisSize,BasisSize);
  end
  
  if (nDim~=BasisSize)
    error('Matrix size (%d) inconsistent with basis size (%d). Please report.',nDim,BasisSize);
  end
  if any(isnan(L))
    error('Liouvillian matrix contains NaN entries! Please report.');
  end
  
  % Rescale by maximum of Hamiltonian superoperator
  if (Opt.Rescale)
    scale = -min(min(imag(L)));
    L = L/scale;
    omega = omega0/scale;
  else
    omega = omega0; % angular frequency
  end
  
  maxDvalLim = 2e3;
  maxDval = max(max(abs(imag(L))));
  logmsg(1,'  size: %dx%d, maxabs: %g',length(L),length(L),full(maxDval));
  if maxDval>maxDvalLim
    error(sprintf('Numerical instability, values in diffusion matrix are too large (%g)!',maxDval));
  end
  
  logmsg(1,'  non-zero elements: %d/%d (%0.2f%%)',nnz(L),length(L).^2,100*nnz(L)/length(L)^2);
  
  %==============================================================
  % Computation of the spectral function
  %==============================================================
  logmsg(1,'Computing spectrum...');
  switch Opt.Solver
    
    case 'L' % Lanczos method
      [alpha,beta,minerr] = chili_lanczos(L,StartingVector,omega,Opt);
      minerr = minerr(end);
      if (minerr<Opt.Threshold)
        thisspec = chili_contfracspec(omega,alpha,beta);
        logmsg(1,'  converged to within %g at iteration %d/%d',...
          Opt.Threshold,numel(alpha),BasisSize);
      else
        thisspec = ones(size(omega));
        logmsg(0,'  Tridiagonalization did not converge to within %g after %d steps!\n  Increase Options.LLKM (current settings [%d,%d,%d,%d])',...
          Opt.Threshold,BasisSize,Opt.LLKM');
      end

    case 'C' % conjugated gradients
      CGshift = 1e-6 + 1e-6i;
      [xx,alpha,beta,err,StepsDone] = chili_conjgrad(L,StartingVector,CGshift);

      logmsg(1,'  step %d/%d: CG converged to within %g',...
        StepsDone,BasisSize,err);

      thisspec = chili_contfracspec(omega,alpha,beta);

    case 'R' % bi-conjugate gradients stabilized
      for iomega = 1:numel(omega)
        u = bicgstab(L+omega(iomega)*speye(size(L)),StartingVector,Opt.Threshold,nDim);
        thisspec(iomega) = real(u'*StartingVector);
      end
      
    case '\' % MATLAB backslash solver for linear system
      I = speye(size(L));
      rho0 = StartingVector;
      for iomega = 1:numel(omega)
        thisspec(iomega) = rho0'*((L+omega(iomega)*I)\rho0);
      end
      thisspec = real(thisspec);
      
    case 'D' % "direct" method by Binsch (eigenbasis)
      L = full(L);
      [U,Lam] = eig(L);
      Lam = diag(Lam);
      rho0 = StartingVector;
      Amplitude = (rho0'*U).'.*(U\rho0);
      thisspec = 0;
      for iPeak = 1:numel(Amplitude)
        thisspec = thisspec + Amplitude(iPeak)./(Lam(iPeak)+omega);
      end

  end

  spec = spec + thisspec*Weights(iOri);
  
end % orientation loop
%==============================================================



%==============================================================
% Phasing
%==============================================================
spec = cos(Exp.mwPhase)*real(spec)+sin(Exp.mwPhase)*imag(spec);
%spec = real(exp(1i*Exp.mwPhase)*spec);
%==============================================================


%==============================================================
% Post-convolution
%==============================================================
if doPostConvolution
  logmsg(1,'Postconvolution...');
  
  % Spin system with shf nuclei only
  pcidx = Opt.PostConvNucs;
  pcSys.g = mean(fullSys.g);
  pcSys.A = mean(fullSys.A(pcidx,:),2);
  if isfield(fullSys,'n')
    pcSys.n = fullSys.n(pcidx);
  end
  pcSys.Nucs = nuclist2string(fullSys.Nucs(pcidx));
  
  % Experimental parameters for isotropic shf spectrum
  if FieldSweep
    pcExp.Range = Exp.Range;
    pcExp.mwFreq = mt2mhz(mean(Exp.Range),pcSys.g)/1e3; % GHz
    Range = Exp.Range;
  else
    pcExp.mwRange = Exp.mwRange;
    pcExp.Field = Exp.Field;
    Range = Exp.mwRange*1e3;
  end
  pcExp.Harmonic = 0;
  pcExp.nPoints = Exp.nPoints;
  
  % Linewidth for shf spectrum
  dx = diff(Range)/(Exp.nPoints-1); % field sweep: mT; freq sweep: MHz
  pcSys.lw = dx/12;
  
  % Simulate isotropic spectrum of shf nuclei
  spec_pc = garlic(pcSys,pcExp);
  spec_pc = spec_pc/sum(spec_pc);
  
  % Convolute SLE spectrum with isotropic spectrum
  spec = conv(spec,spec_pc);
  spec = spec(fix(numel(spec_pc)/2)+(1:Exp.nPoints));
end
%==============================================================




%==============================================================
% Basis set analysis
%==============================================================
Opt.BasisAnalysis = false;
if (Opt.BasisAnalysis)
  logmsg(1,'-------------------------------------------------------------------');
  logmsg(1,'Basis set analysis');
  omega_ = linspace(omega(1),omega(end),12);
  u_sum = 0;
  for iomega = 1:numel(omega_)
    u = bicgstab(L+omega_(iomega)*speye(size(L)),StartingVector,1e-7,180);
    u_sum = u_sum + abs(u)/abs(StartingVector'*u);
  end
  u_sum = u_sum/max(u_sum);

  Thr = [1e-3 1e-4 1e-5 1e-6 1e-8];
  for iThr = 1:numel(Thr)
    inc = (u_sum>Thr(iThr));
    LL = Indices(inc,1);
    Le = max(LL(mod(LL,2)==0));
    Lo = max(LL(mod(LL,2)~=0));
    jK = max(Indices(inc,2));
    K = max(Indices(inc,3));
    M = max(Indices(inc,4));
    fprintf('%0.2e: %3d %3d %3d %3d %3d, size %d\n',Thr(iThr),Le,Lo,jK,K,M,sum(inc));
  end
end
%==============================================================


% Temperature: include Boltzmann equilibrium polarization
%---------------------------------------------------------------
if isfinite(Exp.Temperature)
  if FieldSweep
     DeltaE = planck*Exp.mwFreq*1e9; % joule
  else
     DeltaE = bmagn*mean(Sys.g)*Exp.Field*1e-3; % joule
  end
  e = exp(-DeltaE/boltzm/Exp.Temperature);
  Population = [1 e];
  Population = Population/sum(Population);
  Polarization = Population(1) - Population(2);
  spec = spec*Polarization;
end

% Parallel mode: no intensities
if ParallelMode
  spec = spec*0;
end


%==============================================================
% Convolutional broadening
%---------------------------------------------------------
% Convolution with Gaussian only. Lorentzian broadening is already
% included in the slow-motion simulation via T2.
fwhmG = Sys.lw(1);
if (fwhmG>0) && ConvolutionBroadening
  if FieldSweep
    unitstr = 'mT';
  else
    unitstr = 'MHz';
    fwhmG = fwhmG/1e3; % MHz -> GHz
  end
  dx = xAxis(2) - xAxis(1);
  AlwaysConvolve = true;
  if AlwaysConvolve%(fwhmG/dx>2)
    logmsg(1,'Convoluting with Gaussian (FWHM %g %s)...',fwhmG,unitstr);
    spec = convspec(spec,dx,fwhmG,Exp.ConvHarmonic,1);
    Exp.ConvHarmonic = 0;
  else
    % Skip convolution, since it has no noticeable effect if the linewidth is
    % smaller than about 2*dx.
    Exp.ConvHarmonic = 0;
  end
end

outspec = spec;
%==============================================================


%==============================================================
% Field modulation, or derivatives
%--------------------------------------------------------------
if FieldSweep
  if (Exp.ModAmp>0)
    logmsg(1,'  applying field modulation');
    outspec = fieldmod(xAxis,outspec,Exp.ModAmp,Exp.ModHarmonic);
  else
    if (Exp.DerivHarmonic>0)
      logmsg(1,'  harmonic %d: using differentiation',Exp.DerivHarmonic);
      dx = xAxis(2)-xAxis(1);
      for h = 1:Exp.DerivHarmonic
        dspec = diff(outspec,[],2)/dx;
        outspec = (dspec(:,[1 1:end]) + dspec(:,[1:end end]))/2;
      end
    end
  end
else
  % frequency sweeps: not available
end
%==============================================================



%==============================================================
%  Final processing
%==============================================================

switch (nargout)
case 0,
  cla
  if FieldSweep
    if (xAxis(end)<10000)
      plot(xAxis,outspec);
      xlabel('magnetic field (mT)');
    else
      plot(xAxis/1e3,outspec);
      xlabel('magnetic field (T)');
    end
    axis tight
    ylabel('intensity (arb.u.)');
    title(sprintf('%0.8g GHz, %d points',Exp.mwFreq,numel(xAxis)));
  else
    if (xAxis(end)<1)
      plot(xAxis*1e3,spec);
      xlabel('frequency (MHz)');
    else
      plot(xAxis,spec);
      xlabel('frequency (GHz)');
    end
    axis tight
    ylabel('intensity (arb.u.)');
    title(sprintf('%0.8g mT, %d points',Exp.Field,numel(xAxis)));
  end
case 1,
  varargout = {outspec};
case 2,
  varargout = {xAxis,outspec};
end
%==============================================================


logmsg(1,'-------------------------------------------------------------------');

clear global EasySpinLogLevel

return
%====================================================================
%====================================================================
%====================================================================

%--------------------------------------------------------------------
% Calculates all spin Hamiltonian ISTO coefficients for electron Zeeman,
% hyperfine and nuclear Zeeman terms, for 1 eletron spin-1/2 and up
% to two nuclei.
function Sys = istospinsys(Sys_in,B0,IncludeNZI)

Sys = Sys_in;

% Transformation from molecular frame to diffusion frame
% (DiffFrame contains Euler angles for mol->Diff transformation)
R_M2Diff = erot(Sys.DiffFrame);

% Electron Zeeman
%--------------------------------------------------------------------
for iEl = 1:Sys.nElectrons
  if Sys.fullg
    g = Sys.g(3*(iEl-1)+(1:3),:);
  else
    g = diag(Sys.g(iEl,:));
  end
  if isfield(Sys,'gFrame')
    R_g2M = erot(Sys.gFrame(iEl,:)).';
    g = R_g2M*g*R_g2M.';  % g frame -> molecular frame
  end
  g = R_M2Diff*g*R_M2Diff.';  % molecular frame -> diffusion frame
  
  % Get ISTO components
  [g0,g1,g2] = istocoeff(g);
  if any(abs(g1)>1e-6)
    error('g tensor must be symmetric for this method.');
  end
  
  % Set parameters for chili_liouvmatrix*
  Sys.g_axial(iEl) = g2(1)==0;
  Sys.EZ0(iEl) = bmagn*(B0/1e3)*g0/planck*2*pi; % -> angular frequency
  Sys.EZ2(:,iEl) = bmagn*(B0/1e3)*g2/planck*2*pi; % -> angular frequency
end

% Hyperfine
%--------------------------------------------------------------------
for iEl = 1:Sys.nElectrons
  for iNuc = 1:Sys.nNuclei
    if Sys.fullA
      A = Sys.A(3*(iNuc-1)+(1:3),3*(iEl-1)+(1:3));
    else
      A = diag(Sys.A(iNuc,3*(iEl-1)+(1:3)));
    end
    if isfield(Sys,'AFrame')
      R_M2A = erot(Sys.AFrame(iNuc,3*(iEl-1)+(1:3)));
      R_A2M = R_M2A.';
      A = R_A2M*A*R_A2M.';  % A frame -> molecular frame
    end
    A = R_M2Diff*A*R_M2Diff.';  % molecular frame -> diffusion frame
    
    % Get ISTO components
    [A0,A1,A2] = istocoeff(A);
    if (any(abs(A1)>1e-6))
      error('Hyperfine tensors must be symmetric for this method.');
    end
    
    % Set parameters for chili_liouvmatrix*
    Sys.A_axial(iNuc,iEl) = A2(1)==0;
    Sys.HF0(iNuc,iEl) = A0*1e6*2*pi; % MHz -> angular frequency
    Sys.HF2(:,iNuc,iEl) = A2*1e6*2*pi; % MHz -> angular frequency
  end
end

% Nuclear Zeeman
%--------------------------------------------------------------------
if IncludeNZI
  for iNuc = 1:Sys.nNuclei
    gn0(iNuc) = istocoeff(Sys.gn(iNuc));
    Sys.NZ0(iNuc) = nmagn*(B0/1e3)*gn0(iNuc)/planck*2*pi; % -> angular freq.
  end
else
  Sys.NZ0 = zeros(1,Sys.nNuclei);
end

% Adaption for two nuclei, to feed to chili_liouvmatrix2
if (Sys.nNuclei>=2)
  Sys.Ib = Sys.I(2);
  Sys.NZ0b = Sys.NZ0(2);
  Sys.HF0b = Sys.HF0(2);
  Sys.HF2b = Sys.HF2(:,2);
end

return


%--------------------------------------------------------------------
function Basis = processbasis(Sys,Bas,maxPotentialK)

Basis = Bas;

Basis.evenLmax = Basis.LLKM(1);
Basis.oddLmax = Basis.LLKM(2);
Basis.Kmax = Basis.LLKM(3);
Basis.Mmax = Basis.LLKM(4);

if (Basis.oddLmax>Basis.evenLmax)
  Basis.oddLmax = Basis.evenLmax;
end

% Set jKmin = +1 if tensorial coefficients are all real. This is the
% case when g and A are collinear and the orientation of the diffusion
% tensor is described by the angles (0,beta,0)
if isempty(Basis.jKmin)
  Basis.jKmin = -1;
  if all(isreal(Sys.EZ2))
    if (Sys.nNuclei==0) || (Sys.nNuclei>0) && all(isreal(Sys.HF2))
      Basis.jKmin = +1;
    end
  end
end

if (Sys.nNuclei==0)
  Basis.pImax = 0;
elseif (Sys.nNuclei==1)
  pImax = 2*Sys.I;
  if ~isfield(Basis,'pImax') || isempty(Basis.pImax)
    Basis.pImax = pImax;
  end
  Basis.pImax = min(Basis.pImax,pImax);
else
  pImax = 2*Sys.I;
  if ~isfield(Basis,'pImax') || isempty(Basis.pImax)
    Basis.pImax = pImax;
  end
  Basis.pImax = min(Basis.pImax,pImax);
  Basis.pI1max = Basis.pImax(1);
  Basis.pI2max = Basis.pImax(2);
end

% pSmin
if ~isfield(Basis,'pSmin') || isempty(Basis.pSmin)
  Basis.pSmin = 0;
end

% Use only even K if there is no magnetic or diffusion tilt.
if isempty(Basis.deltaK)
  if (Sys.nNuclei>0)
    noHFtilt = all(Sys.HF2([2 4])==0);
  else
    noHFtilt = true;
  end
  if all(Sys.EZ2([2 4])==0) && noHFtilt
    Basis.deltaK = 2;
  else
    Basis.deltaK = 1;
  end
end

% Use only even L values (oddLmax=0) and no K values (Kmx=0)
% in case of axial magnetic tensors, axial potential, 
% and no magnetic/diffusion tilt
axialSystem = all(Sys.g_axial);
if (Sys.nNuclei>0)
  axialSystem = axialSystem && all(Sys.A_axial(:));
end
if axialSystem && (Basis.deltaK==2) && (maxPotentialK==0)
  Basis.oddLmax = 0;
  Basis.Kmax = 0;
end


Basis.v = [...
  Basis.evenLmax Basis.oddLmax Basis.Kmax Basis.Mmax, ...
  Basis.jKmin Basis.pSmin Basis.deltaK ...
  Basis.MpSymm ...
  Basis.pImax];

return
% 
%========================================================================
function [Dyn,err] = processdynamics(D,FieldSweep)

Dyn = D;
err = '';

% diffusion tensor, correlation time
%------------------------------------------------------------------------
% convert everything (tcorr, logcorr, logDiff) to Diff
if isfield(Dyn,'Diff')
  % Diff given
elseif isfield(Dyn,'logDiff')
  Dyn.Diff = 10.^Dyn.logDiff;
elseif isfield(Dyn,'tcorr')
  Dyn.Diff = 1/6./Dyn.tcorr;
elseif isfield(Dyn,'logtcorr')
  if Dyn.logtcorr>=0, error('Sys.logtcorr must be negative.'); end
  Dyn.Diff = 1/6./10.^Dyn.logtcorr;
else
  err = sprintf('You must specify a rotational correlation time or a diffusion tensor\n(Sys.tcorr, Sys.logtcorr, Sys.Diff or Sys.logDiff).');
  return
end

if any(Dyn.Diff<0)
  error('Negative diffusion rate or correlation times are not possible.');
elseif any(Dyn.Diff>1e12)
  fprintf('Diffusion rate very fast. Simulation might not converge.\n');
elseif any(Dyn.Diff<1e3)
  fprintf('Diffusion rate very slow. Simulation might not converge.\n');
end

% expand to rhombic tensor
switch numel(Dyn.Diff)
  case 1, Dyn.Diff = Dyn.Diff([1 1 1]);
  case 2, Dyn.Diff = Dyn.Diff([1 1 2]);
  case 3, % Diff already rhombic
  otherwise
    err = 'Sys.Diff must have 1, 2 or 3 elements (isotropic, axial, rhombic).';
    return
end

if isfield(Dyn,'lw')
  if numel(Dyn.lw)>1
    if FieldSweep
      LorentzFWHM = Dyn.lw(2)*28 * 1e6; % mT -> MHz -> Hz
    else
      LorentzFWHM = Dyn.lw(2)*1e6; % MHz -> Hz
    end
  else
    LorentzFWHM = 0;
  end
  if (LorentzFWHM~=0)
    % Lorentzian T2 from FWHM in freq domain 1/T2 = pi*FWHM
    Dyn.T2 = 1/LorentzFWHM/pi;
  else
    Dyn.T2 = inf;
  end
end

% Heisenberg exchange
%------------------------------------------------------------------
if ~isfield(Dyn,'Exchange'), Dyn.Exchange = 0; end
Dyn.Exchange = Dyn.Exchange*2*pi*1e6; % MHz -> angular frequency

return
