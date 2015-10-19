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
%     Opt.LLKM            basis size: [evenLmax oddLmax Kmax Mmax]
%     Opt.Verbosity       0: no display, 1: show info
%     Opt.nKnots          number of knots for powder simulation
%     Opt.Symmetry        symmetry to use for powder simulation
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

if Sys.fullg
  mT2MHz = mt2mhz(1,mean(eig(Sys.g)));
else
  mT2MHz = mt2mhz(1,mean(Sys.g));
end
if any(Sys.HStrain) || any(Sys.gStrain) || any(Sys.AStrain) || any(Sys.DStrain)
  error('chili does not support strains (HStrain, gStrain, AStrain, DStrain). Please remove from spin system.');
end

% Convolution with Gaussian only. Lorentzian broadening is 
% included in the slow-motion simulation via T2.
ConvolutionBroadening = any(Sys.lw(1)>0);

% Dynamics
%-------------------------------------------------------------------
if ~isfield(Sys,'lambda'), Sys.lambda = 0; end
if ~isfield(Sys,'Exchange'), Sys.Exchange = 0; end
if ~isfield(Sys,'DiffFrame'), Sys.DiffFrame = [0 0 0]; end

if isfield(Sys,'logtcorr'), Dynamics.logtcorr = Sys.logtcorr; end
if isfield(Sys,'tcorr'), Dynamics.tcorr = Sys.tcorr; end
if isfield(Sys,'logDiff'), Dynamics.logDiff = Sys.logDiff; end
if isfield(Sys,'Diff'), Dynamics.Diff = Sys.Diff; end
if isfield(Sys,'lwpp'), Dynamics.lwpp = Sys.lwpp; end
if isfield(Sys,'lw'), Dynamics.lw = Sys.lw; end
if isfield(Sys,'lambda'), Dynamics.lambda = Sys.lambda; end
if isfield(Sys,'Exchange'), Dynamics.Exchange = Sys.Exchange; end

if isfield(Sys,'psi')
  error('Sys.psi is obsolete. Remove it from your code. See the documentation for details.');
end

% Experimental settings
%-------------------------------------------------------------------
if ~isfield(Exp,'nPoints'), Exp.nPoints = 1024; end
if ~isfield(Exp,'Harmonic'), Exp.Harmonic = []; end
if ~isfield(Exp,'mwPhase'), Exp.mwPhase = 0; end
if ~isfield(Exp,'Temperature'), Exp.Temperature = NaN; end
if ~isfield(Exp,'ModAmp'), Exp.ModAmp = 0; end
if ~isfield(Exp,'Mode'), Exp.Mode = 'perpendicular'; end
if ~isfield(Exp,'Ordering'), Exp.Ordering = []; end
if ~isfield(Exp,'CrystalOrientation'), Exp.CrystalOrientation = []; end

if isfield(Exp,'MOMD')
  error('Exp.MOMD is obsolete. Remove it from your code. See the documentation for details.');
end

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
    logmsg(1,'  partial order (built-in function, lambda = %g)',Exp.Ordering);
  elseif isa(Exp.Ordering,'function_handle')
    UserSuppliedOrderingFcn = true;
    logmsg(1,'  partial order (user-supplied function)');
  else
    error('Exp.Ordering must be a single number or a function handle.');
  end
end

% Determine whether to do a powder simulation
if isempty(Sys.lambda) || all(Sys.lambda==0)
  if isempty(Exp.Ordering) || all(Exp.Ordering==0)
    logmsg(1,'  No ordering potential given, skipping powder simulation.');
    PowderSimulation = false;
  else
  logmsg(1,'  Non-zero ordering potential given, doing powder simulation.');
    PowderSimulation = true;
  end    
else
  if ~isempty(Exp.CrystalOrientation)
    logmsg(1,'  Non-zero ordering potential given, doing single-crystal simulation.');
    PowderSimulation = false;
  else
    logmsg(1,'  Non-zero ordering potential given, doing powder simulation.');
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
if ~isfield(Opt,'Symmetry'), Opt.Symmetry = 'Dinfh'; end
if ~isfield(Opt,'Output'), Opt.Output = 'summed'; end
if ~isfield(Opt,'SymmFrame'), Opt.SymmFrame = []; end
switch Opt.Output
  case 'summed', Opt.SeparateTransitions = 0;
  case 'separate', Opt.SeparateTransitions = 1;
  otherwise, error('Wrong setting in Options.Output.');
end

if Opt.SeparateTransitions
  error('Opt.Output=''separate'' is not supported by chili.');
end

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
  if (Sys.nElectrons==1) && (Sys.S==1/2) && (Sys.nNuclei<=2)
  else
    error('Opt.LiouvMethod=''Freed'' does not work with this spin system.');
  end
end

if generalLiouvillian
  if ~isempty(Sys.lambda) && any(Sys.lambda~=0)
    error('Ordering potential not supported for this spin system.');
  end
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
if ~isfield(Opt,'pImax')
  Opt.pImax = [];
end
Basis.pImax = Opt.pImax;
if ~isfield(Opt,'MeirovitchSymm')
  Opt.MeirovitchSymm = true;
end
Basis.MeirovitchSymm = Opt.MeirovitchSymm;

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

Dynamics.xlk = chili_xlk(Dynamics);
Dynamics.maxL = size(Dynamics.xlk,1)-1;

Basis = processbasis(Sys,Basis,Dynamics);
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
Sys.DirTilt = any(theta~=0);

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

% Set up quantum numbers for basis
%-------------------------------------------------------
logmsg(1,'Setting up orientational basis set...');
logmsg(1,'  Leven max %d, Lodd max %d, Kmax %d, Mmax %d',...
  Basis.LLKM(1),Basis.LLKM(2),Basis.LLKM(3),Basis.LLKM(4));
logmsg(1,'  deltaK %d, jKmin %+d, pSmin %+d, symm %d',...
  Basis.deltaK,Basis.jKmin,Basis.pSmin,Basis.MeirovitchSymm);
if generalLiouvillian
  Basis.List = generatebasis(Basis);
  logmsg(1,'  basis size: %d (%d spatial, %d spin)',size(Basis.List,1)*Sys.nStates^2,size(Basis.List,1),Sys.nStates^2);
else
  [Basis.Size,Basis.SpatialSize,Indices] = chili_basiscount(Basis,Sys);
  %size(Indices), Indices
  logmsg(1,'  basis size: %d (%d spatial, %d spin)',Basis.Size,Basis.SpatialSize,Basis.Size/Basis.SpatialSize);
end

% Preparations
%-----------------------------------------------------------------------
if generalLiouvillian

  % Generate all cartesian spin operators
  for iSpin = 1:numel(Sys.Spins)
    SpinOps{iSpin,1} = sop(Sys.Spins,iSpin,1);
    SpinOps{iSpin,2} = sop(Sys.Spins,iSpin,2);
    SpinOps{iSpin,3} = sop(Sys.Spins,iSpin,3);
  end
  
  [T0,T1,T2,F0,F1,F2] = magint(Sys,SpinOps,CenterField,Opt.IncludeNZI);
  Gamma = rdogamma(Basis.List,Dynamics.Diff,Sys.nStates^2);
  [jjj0,jjj1,jjj2] = jjjsymbol(Basis.LLKM,any(F1(:)));

  % Detection operator
  SxOps = SpinOps{1,1};
  for e = 2:Sys.nElectrons
    SxOps = SxOps + SpinOps{e,1};
  end
  
else
  
  % Pick functions for the calculation of Liouvillian and starting vector
  switch Sys.nNuclei
    case 0
      chili_lm = @chili_lm0;
      chili_sv = @chili_sv0;
    case 1
      chili_lm = @chili_lm1;
      chili_sv = @chili_sv1;
    case 2
      chili_lm = @chili_lm2;
      chili_sv = @chili_sv2;
    otherwise
      error('The chosen method cannot handle %d nuclei.',Sys.nNuclei);
  end
  
end


% Loop over all orientations
%=====================================================================
for iOri = 1:nOrientations
  
  % Set up orientation
  %-------------------------------------------------------
  logmsg(2,'orientation %d of %d: phi = %g°, theta = %g° (weight %g)',...
    iOri,nOrientations,phi(iOri)*180/pi,theta(iOri)*180/pi,Weights(iOri));

  if generalLiouvillian
    D1 = wignerd(1,[phi,theta,0]);
    D2 = wignerd(2,[phi(iOri) theta(iOri) 0]);
    [Q0,Q1,Q2] = rbos(D1,D2,T0,T1,T2,F0,F1,F2);
  else
    Sys.d2psi = wignerd(2,[phi(iOri) theta(iOri) 0]);
  end
  
  % Starting vector
  %-------------------------------------------------------
  logmsg(1,'Computing starting vector(s)...');
  if generalLiouvillian
    StartingVector = startvec(Basis.List,SxOps);
  else
    StartingVector = chili_sv(Sys,Basis,Dynamics,Opt);
  end
  
  BasisSize = size(StartingVector,1);
  nVectors = size(StartingVector,2);
  logmsg(1,'  vector(s) size: %dx%d',BasisSize,nVectors);
  logmsg(1,'  non-zero elements: %d/%d (%0.2f%%)',...
    nnz(StartingVector),numel(StartingVector),100*nnz(StartingVector)/BasisSize);
  logmsg(1,'  maxabs %g, norm %g',full(max(abs(StartingVector))),norm(StartingVector));

  % Liouville matrix
  %-------------------------------------------------------
  logmsg(1,'Computing Liouville matrix...');

  if generalLiouvillian
    H = liouvhamiltonian(Basis.List,Q0,Q1,Q2,jjj0,jjj1,jjj2);
    L = -2i*pi*H + Gamma;
  else
    [r,c,Vals,nDim,nElm] = chili_lm(Sys,Basis.v,Dynamics,Opt.Allocation);
    idx = 1:nElm;
    r = r(idx) + 1;
    c = c(idx) + 1;
    Vals = Vals(idx); % Hz
    if (nDim~=BasisSize)
      error('Matrix size (%d) inconsistent with basis size (%d). Please report.',nDim,BasisSize);
    end
    if any(isnan(Vals))
      error('Liouville matrix contains %d NaN entries!',sum(isnan(Vals)));
    end
    L = sparse(r,c,Vals,BasisSize,BasisSize);
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
  logmsg(1,'  solver: %s',SolverString);
  switch Opt.Solver
    case 'L' % Lanczos method by Jack Freed
      for iVec = 1:nVectors
        [alpha,beta,minerr] = chili_lanczos(L,StartingVector(:,iVec),omega,Opt);
        minerr = minerr(end);
        if (minerr<Opt.Threshold)
          thisspec(iVec,:) = chili_contfracspec(omega,alpha,beta);
          logmsg(1,'  vector %d: converged to within %g at iteration %d/%d',...
            iVec,Opt.Threshold,numel(alpha),BasisSize);
        else
          thisspec = ones(size(omega));
          logmsg(0,'  Tridiagonalization did not converge to within %g after %d steps!\n  Increase Options.LLKM (current settings [%d,%d,%d,%d])',...
            Opt.Threshold,BasisSize,Opt.LLKM');
        end
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
      for iVec = 1:nVectors
        rho0 = StartingVector(:,iVec);
        for iomega = 1:numel(omega)
          thisspec(iVec,iomega) = rho0'*((L+omega(iomega)*I)\rho0);
        end
      end
      thisspec = real(thisspec);
      
    case 'D' %"direct" method by Binsch (eigenbasis)
      L = full(L);
      [U,Lambda] = eig(L);
      Lambda = diag(Lambda);
      thisspec = zeros(nVectors,numel(omega));      
      for iVec = 1:nVectors
        rho0 = StartingVector(:,iVec);
        Amplitude = (rho0'*U).'.*(U\rho0);
        thisspec_ = 0;
        for iPeak = 1:numel(Amplitude)
          thisspec_ = thisspec_ + Amplitude(iPeak)./(Lambda(iPeak)+omega);
        end
        thisspec(iVec,:) = thisspec_;
      end

  end

  for iTrans = 1:nVectors
    spec{iTrans}(iOri,:) = thisspec(iTrans,:);
  end

end % orientation loop
%==============================================================



%==============================================================
% Accumulation of spectra
%==============================================================
% Accumulation
logmsg(1,'Spectra accumulation (%d spectra)',nOrientations);
totalspec = zeros(nVectors,Exp.nPoints);
for iTrans = 1:nVectors
  for iOri = 1:nOrientations
    totalspec(iTrans,:) = totalspec(iTrans,:) + ...
      spec{iTrans}(iOri,:)*Weights(iOri);
  end
end
spec = totalspec;
%==============================================================


%==============================================================
% Phasing
%==============================================================
spec = cos(Exp.mwPhase)*real(spec)+sin(Exp.mwPhase)*imag(spec);
%spec = real(exp(1i*Exp.mwPhase)*spec);
%==============================================================


%==============================================================
% Postconvolution
%==============================================================
if (numel(Sys.I)>2) && ~generalLiouvillian
  logmsg(1,'Postconvolution...');
  shfSys.g = mean(Sys.g);
  shfSys.A = mean(Sys.A(2:end,:),2);
  if isfield(Sys,'n')
    shfSys.n = Sys.n(2:end);
  end
  shfSys.Nucs = nuclist2string(Sys.Nucs(2:end));
  if FieldSweep
    shfExp.Range = Exp.Range;
    shfExp.mwFreq = mt2mhz(mean(Exp.Range),shfSys.g)/1e3; % GHz
    dx = diff(Exp.Range)/(Exp.nPoints-1);
    shfSys.lw = dx/12;
  else
    shfExp.mwRange = Exp.mwRange;
    shfExp.Field = Exp.Field;
    dx = diff(Exp.mwRange*1e3)/(Exp.nPoints-1); % MHz
    shfSys.lw = dx/12;
  end
  shfExp.Harmonic = 0;
  shfExp.nPoints = Exp.nPoints;
  spec_shf = garlic(shfSys,shfExp);
  spec_shf = spec_shf/sum(spec_shf);
  spec = conv(spec,spec_shf);
  spec = spec(fix(numel(spec_shf)/2)+(1:Exp.nPoints));
end
%==============================================================




%==============================================================
% Basis set analysis
%==============================================================
Opt.BasisAnalysis = 0;
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
if Sys.fullg
  g = Sys.g;
else
  g = diag(Sys.g);
end
if isfield(Sys,'gFrame')
  R_g2M = erot(Sys.gFrame);
  g = R_g2M*g*R_g2M.';  % g frame -> diffusion frame
end
g = R_M2Diff*g*R_M2Diff.';  % molecular frame -> diffusion frame

% Get ISTO components
[g0,g1,g2] = istocoeff(g);
if any(abs(g1)>1e-6)
  error('g tensor must be symmetric for this method.');
end

% Set parameters for chili_liouvmatrix*
Sys.g_axial = g2(1)==0;
Sys.EZ0 = bmagn*B0/1e3*g0/planck*2*pi; % -> angular frequency
Sys.EZ2 = bmagn*B0/1e3*g2/planck*2*pi; % -> angular frequency

% Hyperfine
%--------------------------------------------------------------------
for iNuc = 1:Sys.nNuclei
  if Sys.fullA
    A = Sys.A(3*(iNuc-1)+(1:3),:);
  else
    A = diag(Sys.A(iNuc,:));
  end
  if isfield(Sys,'AFrame')
    R_A2M = erot(Sys.AFrame(iNuc,:));
    A = R_A2M*A*R_A2M.';  % A frame -> molecular frame
  end
  A = R_M2Diff*A*R_M2Diff.';  % molecular frame -> diffusion frame
  
  % Get ISTO components
  [A0,A1,A2] = istocoeff(A);
  if (any(abs(A1)>1e-6))
    error('Hyperfine tensors must be symmetric for this method.');
  end
  
  % Set parameters for chili_liouvmatrix*
  Sys.A_axial(iNuc) = A2(1)==0;
  Sys.HF0(iNuc) = A0*1e6*2*pi; % MHz -> angular frequency
  Sys.HF2(:,iNuc) = A2*1e6*2*pi; % MHz -> angular frequency
end

% Nuclear Zeeman
%--------------------------------------------------------------------
for iNuc = 1:Sys.nNuclei
  gn0(iNuc) = istocoeff(Sys.gn(iNuc));
  Sys.NZ0(iNuc) = nmagn*B0/1e3*gn0(iNuc)/planck*2*pi; % -> angular freq.
end
if ~IncludeNZI
  Sys.NZ0 = Sys.NZ0*0;
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
function Basis = processbasis(Sys,Bas,Dyn)

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
axialSystem = Sys.g_axial;
if (Sys.nNuclei>0)
  axialSystem = axialSystem && all(Sys.A_axial);
end
if axialSystem && (Basis.deltaK==2) && (max(Dyn.PotentialK)==0)
  Basis.oddLmax = 0;
  Basis.Kmax = 0;
end


Basis.v = [...
  Basis.evenLmax Basis.oddLmax Basis.Kmax Basis.Mmax, ...
  Basis.jKmin Basis.pSmin Basis.deltaK ...
  Basis.MeirovitchSymm ...
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
if ~isfield(Dyn,'Exchange'), Dyn.Exchange=0; end
Dyn.Exchange = Dyn.Exchange*2*pi*1e6; % MHz -> angular frequency

% Ordering potential
%------------------------------------------------------------------
if ~isfield(Dyn,'lambda'), Dyn.lambda = [0 0 0 0 0]; end
if numel(Dyn.lambda)<5, Dyn.lambda(5) = 0; end
if numel(Dyn.lambda)>5, err = 'Too many potential coefficients!'; return; end

% Process lambda list
Dyn.PotentialL = [2 2 4 4 4];
Dyn.PotentialK = [0 2 0 2 4];
idx = Dyn.lambda~=0;
Dyn.maxL = max(Dyn.PotentialL(idx));
Dyn.maxK = max(Dyn.PotentialK(idx));

return
