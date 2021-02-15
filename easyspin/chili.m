% chili    Simulation of cw EPR spectra in the slow motional regime
%
%   chili(Sys,Exp,Opt)
%   spc = chili(...)
%   [B,spc] = chili(...)
%   [nu,spc] = chili(...)
%
%   Computes a slow-motion cw EPR spectrum.
%
%   Sys: spin system structure
%
%     Sys.tcorr       rotational correlation time (in seconds)
%     Sys.logtcorr    log10 of rotational correlation time (in seconds)
%     Sys.Diff        diffusion rate (s^-1)
%     Sys.logDiff     log10 of diffusion rate (s^-1)
%
%         All fields can have 1 (isotropic), 2 (axial) or 3 (rhombic) elements.
%         Precedence: logtcorr > tcorr > logDiff > Diff.
%
%     Sys.DiffFrame   Euler angles describing the orientation of the
%                     diffusion tensor in the molecular frame (default [0 0 0])
%     Sys.lw          vector with FWHM residual broadenings
%                     1 element:  GaussianFWHM
%                     2 elements: [GaussianFWHM LorentzianFWHM]
%                     units: mT for field sweeps, MHz for frequency sweeps
%     Sys.lwpp        peak-to-peak line widths, same format as Sys.lw
%     Sys.Exchange    spin exchange rate (microsecond^-1)
%     Sys.Potential   orientational potential coefficients
%                       [L1 M1 K1 lambda1; L2 M2 K2 lambda2; ...]
%
%    Exp: experimental parameter settings
%      mwFreq         microwave frequency, in GHz (for field sweeps)
%      Range          sweep range, [sweepmin sweepmax], in mT (for field sweeps)
%      CenterSweep    sweep range, [center sweep], in mT (for field sweeps)
%      Field          static field, in mT (for frequency sweeps)
%      mwRange        sweep range, [sweepmin sweepmax], in GHz (for freq. sweeps)
%      mwCenterSweep  sweep range, [center sweep], in GHz (for freq. sweeps)
%      nPoints        number of points
%      Harmonic       detection harmonic: 0, 1, 2
%      ModAmp         peak-to-peak modulation amplitude, in mT (field sweeps only)
%      mwPhase        detection phase (0 = absorption, pi/2 = dispersion)
%      Temperature    temperature, in K
%
%   Opt: simulation options
%      LLMK           basis set parameters, [evenLmax oddLmax Mmax Kmax]
%      evenK          whether to use only even K values (true/false)
%      highField      whether to use the high-field approximation (true/false)  
%      pImax          maximum nuclear coherence order for basis
%      GridSize       grid size for powder simulation
%      PostConvNucs   nuclei to include perturbationally via post-convolution
%      Verbosity      0: no display, 1: show info
%      Symmetry       symmetry to use for powder simulation
%
%   Output:
%     B               magnetic field axis vector, in mT (for field sweeps)
%     nu              frequency axis vector, in GHz (for frequency sweeps)
%     spc             simulated spectrum, arbitrary units
%
%     If no output arguments are specified, chili plots the simulated spectrum.

function varargout = chili(Sys,Exp,Opt)

if nargin==0, help(mfilename); return; end

% Check expiry date
error(eschecker);

% Check Matlab version
error(chkmlver);

if nargin<2 || nargin>3, error('Wrong number of input arguments!'); end
if nargout<0, error('Not enough output arguments.'); end
if nargout>2, error('Too many output arguments.'); end

if nargin<3, Opt = struct; end

if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 0; % print level
end

global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

%===============================================================================
% Loop over components and isotopologues
%===============================================================================
singleIso = isstruct(Sys) && isfield(Sys,'singleiso') && Sys.singleiso;
if ~singleIso
  logmsg(1,'-- slow motion regime simulation ----------------------------------');
end

FrequencySweep = ~isfield(Exp,'mwFreq') && isfield(Exp,'Field');

if FrequencySweep
  SweepAutoRange = (~isfield(Exp,'mwRange') || isempty(Exp.mwRange)) && ...
    (~isfield(Exp,'mwCenterSweep') || isempty(Exp.mwCenterSweep));
else
  SweepAutoRange = (~isfield(Exp,'Range') || isempty(Exp.Range)) && ...
    (~isfield(Exp,'CenterSweep') || isempty(Exp.CenterSweep));
end
if ~isfield(Opt,'IsoCutoff'), Opt.IsoCutoff = 1e-4; end

if ~isfield(Opt,'Output'), Opt.Output = 'summed'; end
[Output,err] = parseoption(Opt,'Output',{'summed','separate'});
error(err);
summedOutput = Output==1;

if ~isfield(Sys,'singleiso') || ~Sys.singleiso
  
  if ~iscell(Sys), Sys = {Sys}; end
  
  nComponents = numel(Sys);
  logmsg(1,'%d component(s)',nComponents);
  
  % Determine isotopologues for each components
  for c = 1:nComponents
    SysList{c} = isotopologues(Sys{c},Opt.IsoCutoff);
    nIsotopologues(c) = numel(SysList{c});
    logmsg(1,'  component %d: %d isotopologues',c,nIsotopologues(c));
  end
  
  if sum(nIsotopologues)>1 && SweepAutoRange
    if FrequencySweep
      str = 'Exp.mwRange or Exp.mwCenterSweep';
    else
      str = 'Exp.Range or Exp.CenterSweep';
    end
    error('Multiple components: Please specify sweep range manually using %s.',str);
  end
  
  separateSpectra = ~summedOutput && ...
    (nComponents>1 || sum(nIsotopologues)>1);
  if separateSpectra
    spec = [];
    Opt.Output = 'summed'; % summed spectrum for each isotopologue
  else
    spec = 0;
  end
  
  % Loop over all components and isotopologues
  for iComponent = 1:nComponents
    for iIsotopologue = 1:nIsotopologues(iComponent)
      
      % Simulate single-isotopologue spectrum
      Sys_ = SysList{iComponent}(iIsotopologue);
      Sys_.singleiso = true;
      [xAxis,spec_] = chili(Sys_,Exp,Opt);
      
      % Accumulate or append spectra
      if separateSpectra
        spec = [spec; spec_*Sys_.weight];
      else
        spec = spec + spec_*Sys_.weight;
      end
      
    end
  end
  
  % Output and plotting
  switch nargout
    case 0
      plotresults(xAxis,spec,Exp,FrequencySweep);
    case 1
      varargout = {spec};
    case 2
      varargout = {xAxis,spec};
  end
  return
end
%===============================================================================


logmsg(1,'-- component spectrum simulation ----------------------------------');

% Spin system
%-------------------------------------------------------------------------------
if ~isfield(Sys,'Nucs'), Sys.Nucs = ''; end

[Sys,err] = validatespinsys(Sys);
error(err);

% Check for limitations in spin system
if Sys.MO_present, error('chili does not support Sys.Ham* parameters.'); end
if any(Sys.L(:)), error('chili does not support Sys.L.'); end
if isfield(Sys,'nn') && any(Sys.nn(:)~=0)
  error('chili does not support nuclear-nuclear couplings (Sys.nn).');
end
if any(Sys.HStrain(:)) || any(Sys.gStrain(:)) || any(Sys.AStrain(:)) || any(Sys.DStrain(:))
  error('chili does not support strains (HStrain, gStrain, AStrain, DStrain).');
end

% Convolution with Gaussian only. Lorentzian broadening is 
% included in the slow-motion simulation via T2.
ConvolutionBroadening = any(Sys.lw(1)>0);

% Dynamics
%-------------------------------------------------------------------------------
% Add defaults
if ~isfield(Sys,'DiffFrame'), Sys.DiffFrame = [0 0 0]; end
if ~isfield(Sys,'Potential'), Sys.Potential = []; end
if isfield(Sys,'lwpp'), Dynamics.lwpp = Sys.lwpp; end
if isfield(Sys,'lw'), Dynamics.lw = Sys.lw; end


% Heisenberg exchange
%-------------------------------------------------------------------------------
if ~isfield(Sys,'Exchange'), Sys.Exchange = 0; end
Dynamics.Exchange = Sys.Exchange;
Dynamics.Exchange = Dynamics.Exchange*2*pi*1e6; % microseconds^-1 -> rad s^-1

% Check and parse diffusion rate information
%-------------------------------------------------------------------------------
DiffFields = {'logtcorr','tcorr','Diff','logDiff'};
hasDiffFields = isfield(Sys,DiffFields);
if sum(hasDiffFields)>1
  error('Only one of Sys.tcorr, Sys.logtcorr, Sys.Diff, and Sys.logDiff is allowed.');
elseif ~any(hasDiffFields)
  error('One of Sys.tcorr, Sys.logtcorr, Sys.Diff, and Sys.logDiff is required.');
end
for k = 1:numel(DiffFields)
  if hasDiffFields(k)
    Dynamics.(DiffFields{k}) = Sys.(DiffFields{k});
  end
end

% Orientational potential
%-------------------------------------------------------------------------------
% Error on obsolete field Sys.lambda, include explicit upgrade information
if isfield(Sys,'lambda') && ~isempty(Sys.lambda)
  if numel(Sys.lambda)>4
    error('Sys.lambda must be a vector with at most 4 elements.');
  end
  lam = Sys.lambda;
  LMK = [2 0 0; 2 0 2; 4 0 0; 4 0 2];
  str = '    Sys.Potential = [';
  for p = 1:numel(lam)
    str = [str sprintf('%d %d %d %g',LMK(p,1),LMK(p,2),LMK(p,3),lam(p))];
    if p~=numel(lam), str = [str '; ']; end
  end
  str = [str ']'];
  error('\n  Sys.lambda is obsolete.\n  Use the following instead:\n\n%s;    % L M K lambda\n',str);
end

% Extract and organize information about potential
if ~isempty(Sys.Potential)
  if size(Sys.Potential,2)~=4
    error('Sys.Potential needs 4 entries per row (L, M, K, lambda).');
  end
  Potential.L = Sys.Potential(:,1);
  Potential.M = Sys.Potential(:,2);
  Potential.K = Sys.Potential(:,3);
  Potential.lambda = Sys.Potential(:,4);
  rmv = Potential.lambda==0;
  rmv = rmv | (Potential.L==0 & Potential.M==0 & Potential.K==0);
  if any(rmv)
    Potential.L(rmv) = [];
    Potential.M(rmv) = [];
    Potential.K(rmv) = [];
    Potential.lambda(rmv) = [];
  end
else
  Potential.L = [];
  Potential.M = [];
  Potential.K = [];
  Potential.lambda = [];
end
usePotential = ~isempty(Potential.lambda);

% Validate inputs for orientational potential
if usePotential
  if any(Potential.L<0)
    error('L values of potential coefficients must be nonnegative.');
  end
  if any(abs(Potential.K)>Potential.L)
    error('L and K values of potential coefficients do not satisfy |K|<=L.');
  end
  if any(abs(Potential.M)>Potential.L)
    error('L and M values of potential coefficients do not satisfy |M|<=L.');
  end
  if any(Potential.K<0)
    error('Only potential terms with nonnegative values of K are allowed. Terms with negative K required to render the potential real-valued are supplemented automatically.');
  end
  if any(Potential.M(Potential.K==0)<0)
    error('For potential terms with K=0, M must be nonnegative. Terms with negative M required to render the potential real-valued are supplemented automatically.');
  end
  zeroMK = Potential.K==0 & Potential.M==0;
  if any(~isreal(Potential.lambda(zeroMK)))
    error('Potential coefficients for M=K=0 must be real-valued.');
  end
end

% Check for old-style potential (L=2,4; M=0; K=0,2; real-valued lambda)
if usePotential
  oldStylePotential = ...
     all(Potential.L==2 | Potential.L==4) && ...
     all(Potential.M==0) && ...
     all(Potential.K==0 | Potential.K==2) && ...
     all(isreal(Potential.lambda));
end

% Experimental settings
%-------------------------------------------------------------------------------
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

logmsg(1,'Experiment:');

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
  logmsg(1,'  frequency sweep, magnetic field %0.8g T',Exp.Field);
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
  SweepWidth = Exp.CenterSweep(2);
  Exp.Range = Exp.CenterSweep(1) + [-1 1]/2*SweepWidth;
  if any(Exp.Range<0) || diff(Exp.Range)<=0
    error('Invalid sweep range! Check Exp.CenterSweep or Exp.Range.');
  end
else
  CenterFreq = Exp.mwCenterSweep(1);
  SweepWidth = Exp.mwCenterSweep(2);
  Exp.mwRange = Exp.mwCenterSweep(1) + [-1 1]/2*SweepWidth;
  CenterField = Exp.Field;
  if any(Exp.mwRange<0) || diff(Exp.mwRange)<=0
    error('Invalid sweep range! Check Exp.mwCenterSweep or Exp.mwRange.');
  end
end

if FieldSweep
  logmsg(1,'  field range (mT): min %g, max %g, center %g, width %g',...
    Exp.Range(1),Exp.Range(2),CenterField,SweepWidth);
else
  logmsg(1,'  frequency range (GHz): min %g, max %g, center %g, width %g',...
    Exp.mwRange(1),Exp.mwRange(2),CenterFreq,SweepWidth);
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
if Exp.ModAmp>0
  if FieldSweep
    logmsg(1,'  field modulation, amplitude %g mT',Exp.ModAmp);
    if Exp.Harmonic<1
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
if ParallelMode
  error('chili does not support parallel-mode spectra.');
end

logmsg(1,'  %d points',Exp.nPoints);

% Complain if fields only valid in pepper() are given
if isfield(Exp,'Orientations')
  warning('Exp.Orientations is obsolete. Use Exp.CrystalOrientations instead.');
end
if isfield(Exp,'CrystalSymmetry')
  warning('Exp.CrystalSymmetry is not used by chili.');
end

% Ordering of director frame (partial ordering)
if ~isempty(Exp.Ordering)
  if isnumeric(Exp.Ordering) && numel(Exp.Ordering)==1 && isreal(Exp.Ordering)
    lam = Exp.Ordering;
    if lam~=0
      Exp.Ordering = @(phi,theta) exp(lam*plegendre(2,0,cos(theta))).*ones(size(phi));
      logmsg(1,'  director ordering: built-in function, coefficient = %g',lam);
    else
      Exp.Ordering = [];
      logmsg(1,'  director ordering: none');
    end
  elseif isa(Exp.Ordering,'function_handle')
    if nargin(Exp.Ordering)<2
      error('The function in Exp.Ordering must accept two inputs.');
    end
    if nargout(Exp.Ordering)<1
      error('The function in Exp.Ordering must provide one output.');
    end
    logmsg(1,'  director ordering: user-supplied function)');
  else
    error('Exp.Ordering must be a single number or a function handle.');
  end
else
  logmsg(1,'  director ordering: none');
end
useDirectorOrdering = ~isempty(Exp.Ordering);

% Determine whether to do a powder simulation
% (without potential, no powder sim is necessary - it's identical to a
% single-orientation sim)
PowderSimulation = isempty(Exp.CrystalOrientation) && usePotential;

% Options
%-------------------------------------------------------------------------------
if isempty(Opt), Opt = struct; end

% Documented
if ~isfield(Opt,'LLMK')
  if usePotential
    error('Sys.Potential is given. Please provide basis set information (Opt.LLMK etc).');
  else
    Opt.LLMK = [14 7 2 6];
  end
end
if ~isfield(Opt,'highField'), Opt.highField = false; end
if ~isfield(Opt,'pImax'), Opt.pImax = []; end
if ~isfield(Opt,'pImaxall'), Opt.pImaxall = []; end
if ~isfield(Opt,'GridSize'), Opt.GridSize = [19 0]; end
if ~isfield(Opt,'LiouvMethod'), Opt.LiouvMethod = ''; end
if ~isfield(Opt,'FieldSweepMethod'), Opt.FieldSweepMethod = []; end
if ~isfield(Opt,'PostConvNucs'), Opt.PostConvNucs = ''; end
if ~isfield(Opt,'Solver'), Opt.Solver = ''; end
if ~isfield(Opt,'GridSymmetry'), Opt.GridSymmetry = 'Dinfh'; end
% Opt.Verbosity

% Undocumented
if ~isfield(Opt,'Rescale'), Opt.Rescale = true; end
if ~isfield(Opt,'Threshold'), Opt.Threshold = 1e-6; end
if ~isfield(Opt,'Lentz'), Opt.Lentz = true; end
if ~isfield(Opt,'IncludeNZI'), Opt.IncludeNZI = true; end
if ~isfield(Opt,'pqOrder'), Opt.pqOrder = false; end
if ~isfield(Opt,'GridFrame'), Opt.GridFrame = []; end
if ~isfield(Opt,'Diagnostics'), Opt.Diagnostics = ''; end
if ~isfield(Opt,'useLMKbasis'), Opt.useLMKbasis = false; end
if ~isfield(Opt,'useStartvecSelectionRules'), Opt.useStartvecSelectionRules = true; end
if ~isfield(Opt,'PeqTol'), Opt.PeqTol = []; end
if ~isfield(Opt,'jKmin'), Opt.jKmin = []; end
if ~isfield(Opt,'evenK'), Opt.evenK = []; end

if isfield(Opt,'deltaK')
  error('Opt.deltaK is obsolete. Use Opt.evenK instead.');
end

if ~ischar(Opt.Diagnostics) && ~isempty(Opt.Diagnostics) && ~isvarname(Opt.Diagnostics)
  error('If given, Opt.Diagnosics must be a valid Matlab variable name.');
end
saveDiagnostics = ~isempty(Opt.Diagnostics);

% Determine default method for constructing Liouvillian
if ~isfield(Opt,'LiouvMethod') || isempty(Opt.LiouvMethod)
  if (Sys.nElectrons==1) && (Sys.S==1/2) && (Sys.nNuclei<=2) && ...
      (~usePotential || oldStylePotential)
    Opt.LiouvMethod = 'fast';
  else
    Opt.LiouvMethod = 'general';
  end
end

[LiouvMethod,err] = parseoption(Opt,'LiouvMethod',{'fast','general'});
error(err);
generalLiouvillian = LiouvMethod==2;

if generalLiouvillian
  if any(Sys.Exchange~=0)
    error('Opt.LiouvMethod=''general'' does not support spin exchange (Sys.Exchange).');
  end
else
  if Sys.nElectrons>1 || Sys.S~=1/2 || Sys.nNuclei>2
    error('Opt.LiouvMethod=''fast'' does not work with this spin system.');
  end
  if usePotential
    if ~oldStylePotential
      error('Opt.LiouvMethod=''fast'' does not work with this orientational potential.');
    end
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
    error('Cannot have more than two nuclei for the Stochastic Liouville equation with this Opt.LiouvMethod.');
  end
  fullSys = Sys;
  Sys = nucspinrmv(Sys,Opt.PostConvNucs);
  Sys.processed = false;
  [Sys,err] = validatespinsys(Sys);
  error(err);
end

if any(Sys.n~=1)
  error('chili cannot handle systems with nuclei with Sys.n > 1 only if these nuclei are treated using post-convolution (Opt.PostConvNucs).');
end

if numel(Opt.GridSize)<1, Opt.GridSize(1) = 19; end
if numel(Opt.GridSize)<2, Opt.GridSize(2) = 0; end
if Opt.GridSize(2)~=0
  error('chili cannot interpolate orientations. Set Opt.GridSize(2) to zero.');
end

% Basis settings
%-------------------------------------------------------------------------------
if isfield(Opt,'LLKM') % error if pre-6.0 field name is given
  LLKM = Opt.LLKM;
  error(['Opt.LLKM is obsolete. Use Opt.LLMK instead. Note that M and K are swapped.\n' ...
    'To convert Opt.LLKM to Opt.LLMK, swap the third and the fourth number.\n' ...
    'In the current case, use\n' ...
    '  Opt.LLMK = [%d %d %d %d];'],LLKM(1),LLKM(2),LLKM(4),LLKM(3));
end
if ~isnumeric(Opt.LLMK) || numel(Opt.LLMK)~=4 || any(Opt.LLMK<0) || any(mod(Opt.LLMK,1))
  error('Opt.LLMK must be a 4-element array with non-negative integers.');
end
if mod(Opt.LLMK(1),2)
  error('Opt.LLMK(1) must be an even number.');
end
if Opt.LLMK(2)>0 && ~mod(Opt.LLMK(2),2)
  error('Opt.LLMK(2) must be an odd integer.');
end
maxL = max(Opt.LLMK(1:2));
if any(Opt.LLMK(3)>maxL)
  error('Opt.LLMK(3)=%d, but must be smaller or equal to the maximum L (%d).',Opt.LLMK(3),maxL);
end
if any(Opt.LLMK(4)>maxL)
  error('Opt.LLMK(4)=%d, but must be smaller or equal to the maximum L (%d).',Opt.LLMK(4),maxL);
end
Basis.LLMK = Opt.LLMK;
Basis.jKmin = Opt.jKmin;
Basis.evenK = Opt.evenK;

% high-field approximation
if isfield(Opt,'pSmin')
  error('Opt.pSmin is not supported - use Opt.highField instead.');
end
if Opt.highField
  Basis.pSmin = +1;
else
  Basis.pSmin = -1;
end

% pImax (maximum nuclear coherence order)
if ~isempty(Opt.pImax)
  if numel(Opt.pImax)~=Sys.nNuclei && numel(Opt.pImax)~=1
    error('Opt.pImax must contain one entry for every nucleus or just a single number.');
  end
  if any(Opt.pImax<0)
    error('Every element in Opt.pImax must be 0 or larger.');
  end
end
Basis.pImax = Opt.pImax;
if ~isempty(Opt.pImaxall)
  if numel(Opt.pImaxall)~=1
    error('Opt.pImaxall must be a single non-negative number.');
  end
end
Basis.pImaxall = Opt.pImaxall;

% M-pS-pI symmetry (see Meirovitch J.Chem.Phys. 77 3915 (1982), eq. (A47))
if ~isfield(Opt,'MpSymm')
  Opt.MpSymm = false;
end
Basis.MpSymm = Opt.MpSymm;


% Field sweep and linear solver
%-------------------------------------------------------------------------------
logmsg(1,'Methods:');

% Field sweep method: set default
if FieldSweep
  if isempty(Opt.FieldSweepMethod)
    Opt.FieldSweepMethod = 'approxlin';
  end
  [~,err] = parseoption(Opt,'FieldSweepMethod',{'explicit','approxinv','approxlin'});
  error(err);
end
explicitFieldSweep = FieldSweep && strcmp(Opt.FieldSweepMethod,'explicit');
logmsg(1,'  field sweep method: %s',Opt.FieldSweepMethod);

% Set solver if not given
if isempty(Opt.Solver)
  if strcmp(Opt.FieldSweepMethod,'explicit')
    Opt.Solver = '\';
  else
    Opt.Solver = 'L';
  end
end
useLanczosSolver = Opt.Solver=='L';

if explicitFieldSweep && Opt.Solver~='\'
  error('For an explicit field sweep, use Opt.Solver=''\''.');
end

switch Opt.Solver
  case 'L'
    if Opt.Lentz==1 % Lentz method
      SolverString = 'Lanczos tridiagonalization, left-to-right continued fraction evaluation';
    else
      SolverString = 'Lanczos tridiagonalization, right-to-left continued fraction evaluation';
    end
  case '\'
    SolverString = 'backslash';
  case 'E'
    SolverString = 'eigenvalue method, sum of Lorentzians';
  case 'C'
    SolverString = 'conjugate gradients tridiagonalization, right-to-left continued fraction evaluation';
  case 'B'
    SolverString = 'biconjugate gradients, stabilized';
  otherwise
    error('Unknown method in Options.Solver. Must be ''L'', ''\'', ''E'', ''C'', or ''B''.');
end
logmsg(1,'  solver: %s = %s',Opt.Solver,SolverString);


% Set allocation block size, used in chili_lm
%-------------------------------------------------------------------------------
if ~generalLiouvillian
  blockSize = 1e6;
  minBlockSize = 1e3;
  if ~isfield(Opt,'AllocationBlockSize')
    Opt.AllocationBlockSize = blockSize;
  end
  Opt.AllocationBlockSize = Opt.AllocationBlockSize(1);
  if Opt.AllocationBlockSize<minBlockSize
    error('Opt.AllocationBlockSize = %d is too small. Increase its value to at least %d.',Opt.AllocationBlockSize,minBlockSize);
  end
  logmsg(2,'  allocating memory in blocks of %d non-zero elements',Opt.AllocationBlockSize);
end

% Precalculate spin operator matrices
%-------------------------------------------------------------------------------
if generalLiouvillian
  logmsg(1,'  construction of Liouvillian: general code');
  
  % calculate spin operators
  for iSpin = 1:numel(Sys.Spins)
    SpinOps{iSpin,1} = sop(Sys.Spins,[iSpin,1],'sparse'); % Sx
    SpinOps{iSpin,2} = sop(Sys.Spins,[iSpin,2],'sparse'); % Sy
    SpinOps{iSpin,3} = sop(Sys.Spins,[iSpin,3],'sparse'); % Sz
  end

  SdetOp = sparse(0);
  for e = 1:Sys.nElectrons
    SdetOp = SdetOp + SpinOps{e,1}; % Sx
  end
  
else
  logmsg(1,'  Liouvillian construction: fast S=1/2 code');
  % no need to calculate spin operators for the fast code
  SpinOps = [];
end

% Calculate ISTOs and symmetry properties
%-------------------------------------------------------------------------------
includeNQI = false;
[T,F,Sys,Symmetry,isFieldDep] = magint(Sys,SpinOps,CenterField,...
                                       Opt.IncludeNZI,includeNQI,...
                                       explicitFieldSweep);

if saveDiagnostics
  diagnostics.T = T;
  diagnostics.F = F;
end
                                     
noAnisotropiesPresent = all(F.F1(:)==0) && all(F.F2(:)==0);
if noAnisotropiesPresent
  error('This is an isotropic spin system. chili cannot calculate a slow-motion spectrum.');
end

% Process diffusion tensor and linewidth
%-------------------------------------------------------------------------------
[Dynamics,err] = processdynamics(Dynamics,FieldSweep);
error(err);


% Basis
%-------------------------------------------------------------------------------
Basis = processbasis(Basis,max(Potential.K),Sys.I,Symmetry);


% Set up g values
%-------------------------------------------------------------------------------
% Set reference g value (for frequency-to-field conversion and Boltzmann)
if Sys.fullg
  idx = 1:3;
  for iElectron = 1:Sys.nElectrons
    mean_g(iElectron) = mean(eig(Sys.g(idx,:)));
    idx = idx + 3;
  end
  gavg = mean(mean_g);
else
  gavg = mean(mean(Sys.g));
end

% Set field and frequency axes and values
%-------------------------------------------------------------------------------
% Determine B0 and nu based on sweep domain and method
if FieldSweep
  B_ = linspace(Exp.Range(1),Exp.Range(2),Exp.nPoints);
  switch Opt.FieldSweepMethod
    case 'explicit'
      B0 = B_;
      nu = Exp.mwFreq;
    case 'approxinv'
      B0 = CenterField;
      nu = Exp.mwFreq*B0./B_;
    case 'approxlin'
      B0 = CenterField;
      nu = Exp.mwFreq - mt2mhz(B_-B0,gavg)/1e3; % GHz
    otherwise
      error('Unknown setting ''%s'' in Opt.FieldSweepMethod.',Opt.FieldSweepMethod);
  end
else % frequency sweep
  B0 = CenterField;
  nu = linspace(Exp.mwRange(1),Exp.mwRange(2),Exp.nPoints);
end

% Adjust units
B0 = B0/1e3; % mT -> T
omega0 = complex(2*pi*nu*1e9,1/Dynamics.T2); % GHz -> rad s^-1

% Set x axis
if FieldSweep
  xAxis = linspace(Exp.Range(1),Exp.Range(2),Exp.nPoints);  % field axis, mT
else
  xAxis = nu; % frequency axis, GHz
end


% Set up list of orientations
%===============================================================================
if PowderSimulation
  if Opt.GridSize(1)==1
    phi = 0;
    theta = 0;
    GridWeights = 4*pi;
  else
    grid = sphgrid(Opt.GridSymmetry,Opt.GridSize(1));
    Vecs = grid.vecs;
    GridWeights = grid.weights;
    % Transform vector to reference frame representation and convert to polar angles.
    if isempty(Opt.GridFrame)
      [phi,theta] = vec2ang(Vecs);
    else
      [phi,theta] = vec2ang(Opt.GridFrame*Vecs);
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
if useDirectorOrdering
  orifun = foldoridist(Exp.Ordering,Opt.GridSymmetry);
  OrderingWeights = orifun(phi,theta);
  if any(OrderingWeights<0), error('User-supplied orientation distribution gives negative values!'); end
  if all(OrderingWeights==0), error('User-supplied orientation distribution is all-zero.'); end
  logmsg(2,'  orientational potential');
else
  OrderingWeights = ones(1,nOrientations);
end

Weights = GridWeights.*OrderingWeights;
Weights = 4*pi*Weights/sum(Weights);


% Basis set preparations
%-------------------------------------------------------------------------------
logmsg(1,'Basis set:');
logmsg(1,'  truncation settings:');
logmsg(1,'    orientational: Leven max %d, Lodd max %d, Mmax %d, Kmax %d, evenK %d, jKmin %+d',...
  Basis.LLMK(1),Basis.LLMK(2),Basis.LLMK(3),Basis.LLMK(4),Basis.evenK,Basis.jKmin);
logmsg(1,'    spin: pSmin %+d, pImax %s, pImaxall %d',Basis.pSmin,num2str(Basis.pImax),Basis.pImaxall);
logmsg(1,'    M-p symmetry: %d',Basis.MpSymm);

if generalLiouvillian
  
  % Set up basis
  if Opt.useLMKbasis
    Basis = generateoribasis(Basis,'LMK');
  else
    Basis = generateoribasis(Basis,'LjKKM');
  end
  nOriBasis = numel(Basis.L);
  nSpinBasis = Sys.nStates^2;
  logmsg(1,'  orientational basis: %d functions',nOriBasis);
  
  % Get (p,q) quantum numbers for transitions, and index vector for reordering
  % basis states from m1-m2 order (standard) to p-q order (used in fast code)
  [idxpq,pq] = pqorder(Sys.Spins);
  % Reorder detection operator if needed
  SdetOp = SdetOp(:);
  if Opt.pqOrder
    SdetOp = SdetOp(idxpq);
  end
  
  % Removing unwanted spin functions (in mm ordering)
  keep = true;
  % (1) keep only transitions with pS>=pSmin, for each electron
  for ie = 1:Sys.nElectrons
    pS = pq(:,2*ie-1);
    keep = keep & (pS>=Basis.pSmin);
  end
  % (2a) keep only transitions with |pI|<=pImax, for each nucleus
  % (2b) keep only transitions with sum(pI)<=pImaxall
  pIsum = 0;
  for in = 1:Sys.nNuclei
    pI = pq(:,2*Sys.nElectrons+2*in-1);
    pIsum = pIsum + pI;
    keep = keep & (abs(pI)<=Basis.pImax(in));
  end
  keep = keep & (abs(pIsum)<=Basis.pImaxall);
  if Opt.pqOrder
    keep = keep(idxpq);
  end
  logmsg(1,'  spin basis: %d functions (%0.1f%% of %d)',...
    sum(keep),sum(keep)/nSpinBasis*100,nSpinBasis);
  keep = repmat(keep,nOriBasis,1);
  
  % Apply M=p-1 symmetry (Meirovitch J.Chem.Phys. 77(7), 3915-3938, Eq. (A47))
  if Opt.MpSymm
    M = Basis.M;
    psum = sum(pq(:,1:2:end),2);
    keep_Mp = bsxfun(@minus,psum,M.')==1; % keep only basis states with pS+pI-M == 1
    keep = keep & keep_Mp(:);
    logmsg(1,'  applying M-p symmetry: keeping %d of %d functions',sum(keep),numel(keep));
  end

  logmsg(1,'  total basis size: %d functions',sum(keep));
  
else
  
  Basis = chili_basisbuild(Basis,Sys);
  logmsg(1,'  total basis size: %d functions',numel(Basis.L));
  
end
if saveDiagnostics
  diagnostics.basis = Basis;
end

% Precalculate 3j symbols
%-------------------------------------------------------------------------------
if generalLiouvillian
  
  logmsg(1,'Precalculating 3j symbols up to Leven=%d and Lodd=%d...',...
    Basis.LLMK(1),Basis.LLMK(2));
  computeRankOne = any(F.F1(:)~=0);
  [jjj0,jjj1,jjj2] = jjjsymbol(Basis.LLMK(1),Basis.LLMK(2),computeRankOne);
  
end

% Calculate diffusion operator matrix
%-------------------------------------------------------------------------------
% Pre-calculate diffusion operator Wigner expansion coefficient
% (needed for both Opt.LiouvMethod='fast' and 'general')
if usePotential
  logmsg(1,'Calculating Wigner expansion coefficients for potential-dependent part of diffusion matrix');
  XLMK = chili_xlmk(Potential,Dynamics.R);
else
  XLMK = {};
end

if generalLiouvillian

  logmsg(1,'Calculating the diffusion matrix');
  
  % Calculate diffusion superoperator in spatial basis
  Gamma = diffsuperop(Basis,Dynamics.R,XLMK,Potential);
  % Expand to full product basis
  Gamma = kroneye(Gamma,Sys.nStates^2);
  Gamma = Gamma(keep,keep); % prune
  
  maxerr = @(x) full(max(abs(x(:))));
  imagerr = @(x) maxerr(imag(x))/maxerr(real(x));
  logmsg(1,'  imag/real = %f',imagerr(Gamma));
  
else
  
  % Gamma is calculated simultaneously with H
  
end


% Starting vector
%-------------------------------------------------------------------------------
logmsg(1,'Calculating starting vector...');
if generalLiouvillian
  if ~isfield(Opt,'StartVec') || isempty(Opt.StartVec)
    % Calculate sqrt(Peq) vector
    Opt_.useSelectionRules = Opt.useStartvecSelectionRules;
    Opt_.PeqTolerances = Opt.PeqTol;
    [sqrtPeq,nInt] = chili_eqpopvec(Basis,Potential,Opt_);
    % Set up in full product basis, then prune
    StartVector = kron(sqrtPeq,SdetOp(:)/norm(SdetOp(:)));
    StartVector = StartVector(keep);
    normPeqVec = norm(sqrtPeq)^2;
    logmsg(1,'  norm of Peq vector: %g',normPeqVec);
    if normPeqVec<0.99
      fprintf('The norm of the equilibrium population vector in this truncated basis is %g. It should be close to 1. The basis might be too small.',normPeqVec);
    end
  else
    logmsg(1,'  using provided vector');
    if numel(Opt.StartVec)==sum(keep)
      StartVector = Opt.StartVec;
      nInt = [];
    else
      error('Opt.StartVec must have %d elements.',sum(keep));
    end
  end
  StartVector = StartVector/norm(StartVector);
else
  if usePotential
    % Organize potential as expected by chili_startingvector
    LMK = [2 0 0; 2 0 2; 4 0 0; 4 0 2]; % standard terms and order for fast code
    PotLMK = [Potential.L Potential.M Potential.K];
    lambda_shortlist = [0; 0; 0; 0];
    for p = 1:4
      [found,idx] = ismember(LMK(p,:),PotLMK,'rows');
      if found, lambda_shortlist(p) = Potential.lambda(idx); end
    end
  else
    lambda_shortlist = [];
  end
  StartVector = chili_startingvector(Basis,lambda_shortlist);
  nInt = [];
end

if saveDiagnostics
  diagnostics.sv = StartVector;
end

BasisSize = size(StartVector,1);
logmsg(1,'  vector size: %dx1',BasisSize);
logmsg(1,'  non-zero elements: %d/%d (%0.2f%%)',...
  nnz(StartVector),BasisSize,100*nnz(StartVector)/BasisSize);
logmsg(1,'  maxabs %g, norm %g',full(max(abs(StartVector))),norm(StartVector));
if ~isempty(nInt)
  logmsg(1,'  evaluated integrals: 1D %d, 2D %d, 3D %d',nInt(1),nInt(2),nInt(3));
end

% Prepare transformation matrix for K-symmetrization
if Opt.useLMKbasis && useLanczosSolver
  TT = ksymmetrizer(Basis); % in spatial basis
  TT = kron(TT,eye(Sys.nStates^2)); % expand to full basis, incl. spin
  TT = TT(keep,keep); % prune
end


% Loop over all orientations
%===============================================================================
spec = 0;
for iOri = 1:nOrientations
  
  logmsg(1,'Orientation %d/%d: phi = %gdeg, theta = %gdeg (weight %g)',...
    iOri,nOrientations,phi(iOri)*180/pi,theta(iOri)*180/pi,Weights(iOri));

  % Liouvillian matrix
  %-----------------------------------------------------------------------------
  logmsg(1,'Computing Liouvillian matrix...');
  if generalLiouvillian
    [Q0B,Q1B,Q2B,Q0G,Q1G,Q2G] = rbos(T,F,[phi(iOri),theta(iOri),0],isFieldDep);
    
    if Opt.pqOrder
      Q0B = Q0B(idxpq,idxpq);
      Q0G = Q0G(idxpq,idxpq);
      for k = 1:numel(Q1B)
        Q1B{k} = Q1B{k}(idxpq,idxpq);
        Q1G{k} = Q1G{k}(idxpq,idxpq);
      end
      for k = 1:numel(Q2B)
        Q2B{k} = Q2B{k}(idxpq,idxpq);
        Q2G{k} = Q2G{k}(idxpq,idxpq);
      end
    end
    if explicitFieldSweep
      if Opt.useLMKbasis
        HB = liouvhamiltonian_LMK(Basis,Q0B,Q1B,Q2B,jjj0,jjj1,jjj2);
        HG = liouvhamiltonian_LMK(Basis,Q0G,Q1G,Q2G,jjj0,jjj1,jjj2);
      else
        HB = liouvhamiltonian(Basis,Q0B,Q1B,Q2B,jjj0,jjj1,jjj2);
        HG = liouvhamiltonian(Basis,Q0G,Q1G,Q2G,jjj0,jjj1,jjj2);
      end
      HB = HB(keep,keep);
      HG = HG(keep,keep);
    else
      Q0 = Q0B+Q0G;
      if any(F.F1(:))
        for i = 1:3
          for j = 1:3
            Q1{i,j} = Q1B{i,j} + Q1G{i,j};
          end
        end
      else
        Q1 = {};
      end
      for i = 1:5
        for j = 1:5
          Q2{i,j} = Q2B{i,j} + Q2G{i,j};
        end
      end
      if Opt.useLMKbasis
        H = liouvhamiltonian_LMK(Basis,Q0,Q1,Q2,jjj0,jjj1,jjj2);
      else
        H = liouvhamiltonian(Basis,Q0,Q1,Q2,jjj0,jjj1,jjj2);
      end
      H = H(keep,keep);
      if saveDiagnostics
        diagnostics.Q0 = Q0;
        diagnostics.Q1 = Q1;
        diagnostics.Q2 = Q2;        
      end
    end
    
  else
    
    Sys.d2psi = wignerd(2,phi(iOri),theta(iOri),0);
    if explicitFieldSweep
      EZ0_ = Sys.EZ0;
      EZ2_ = Sys.EZ2;
      if isfield(Sys,'NZ0')
        for iNuc = 1:numel(Sys.NZ0)
          NZ0_(iNuc) = Sys.NZ0(iNuc);
        end
      end
    end
    
  end
    
  % Loop over all field values
  %-----------------------------------------------------------------------------
  for iB = 1:numel(B0)
    logmsg(1,'Field value %d/%d: %g mT',iB,numel(B0),B0(iB));
    
    if generalLiouvillian
      
      if explicitFieldSweep
        H = B0(iB)*HB + HG;
      end
      L = 2i*pi*H + Gamma;  % Hamiltonian: Hz -> rad s^-1
      nDim = size(L,1);
      
      if nDim~=BasisSize
        error('Matrix size (%d) inconsistent with basis size (%d). Please report.',nDim,BasisSize);
      end
      if any(isnan(L))
        error('Liouvillian matrix contains NaN entries! Please report.');
      end
      
    else
      
      if explicitFieldSweep
        Sys.EZ0 = EZ0_*B0(iB);
        Sys.EZ2 = EZ2_*B0(iB);
        if isfield(Sys,'NZ0')
          for iNuc = 1:numel(Sys.NZ0)
            Sys.NZ0(iNuc) = NZ0_(iNuc)*B0(iB);
          end
        end
      end
      Sys.DirTilt = Basis.DirTilt; % used in chili_lm
      
      % Build xlk array needed by chili_lm (rearranged from XLMK)
      maxL = numel(XLMK)-1; % maxmimum L in XLMK ( = 2*L from potential)
      xlk = [];
      for L_ = 0:maxL
        xlk(L_+1,1:2*L_+1) = XLMK{L_+1}(L_+1,:);
      end
      Dynamics.xlk = xlk;
      Dynamics.maxL = maxL;      
      Dynamics.Diff = Dynamics.R;
      
      % Call mex function to get L matrix elements
      [r,c,Vals,nDim] = chili_lm(Sys,Basis.v,Dynamics,Opt.AllocationBlockSize);
      % (chili_lm constructs r/c/Vals for -1i*H + Gamma)
      Vals = conj(Vals); % make sure L = +1i*H + Gamma (assumes Gamma is real)
      L = sparse(r,c,Vals,BasisSize,BasisSize);
      
      if saveDiagnostics && iOri==1
        % extract H and Gamma from L = 1i*H + Gamma
        % (assumes only that H and Gamma are Hermitian)
        reL = real(L);
        imL = imag(L);
        Gamma = (reL+reL.')/2 + 1i*(imL-imL.')/2;
        H = (imL+imL.')/2 + 1i*(reL.'-reL)/2;
        H = H/(2*pi); % rad s^-1 -> Hz
      end
      
    end
    
    if saveDiagnostics && iOri==1
      diagnostics.iOri = iOri;
      diagnostics.L = L;
      diagnostics.H = H;
      diagnostics.Gamma = Gamma;
    end
    
    % Apply K-symmetrization if needed to obtain complex symmetric L for Lanczos
    % algorithm. L = i*H + Gamma is complex symmetric only if both H and Gamma
    % are real-valued.
    if generalLiouvillian && Opt.useLMKbasis && useLanczosSolver
      isComplexSymmetric = isreal(H);
      %maxerr = @(A)max(abs(A(:)));
      %imagerr = @(A)maxerr(imag(A))/maxerr(real(A));
      %Himag = imagerr(H)
      %Gimag = imagerr(Gamma)
      if ~isComplexSymmetric
        L = TT'*L*TT;
        %ksymmHimag = imagerr(TT'*H*TT);
        %ksymmHimag = imagerr(TT'*Gamma*TT);
        StartVector = TT'*StartVector;
      end
    end
    
    % Rescale for numerical stability
    if Opt.Rescale
      scale = max(abs(L(:)));
      L = L/scale;
      omega = omega0/scale;
    else
      scale = 1;
      omega = omega0;
    end
    
    maxDval = max(abs(real(L(:))));
    logmsg(2,'  size: %dx%d, maxabsreal: %g',length(L),length(L),full(maxDval));
    
    maxDvalLim = 2e3;
    if maxDval>maxDvalLim
      %  error(sprintf('Numerical instability, values in diffusion matrix are too large (%g)!',maxDval));
    end
    
    logmsg(2,'  non-zero elements: %d/%d (%0.2f%%)',nnz(L),length(L).^2,100*nnz(L)/length(L)^2);
    
    
    %===========================================================================
    % Computation of the spectral function
    %===========================================================================
    switch Opt.Solver
      
      case 'L' % Lanczos method
        if generalLiouvillian && usePotential
          maxabs = @(a)max(abs(a(:)));
          isComplexSymmetric = maxabs(L-L.')/maxabs(L) < 1e-10;
          if ~isComplexSymmetric
            error('L is not complex symmetric - cannot use Lanczos method.');
          end
        end
        [thisspec,converged,dspec] = chili_lanczos(L,StartVector,-1i*omega,Opt);
        if converged
          logmsg(2,'  converged to within %g at iteration %d/%d',...
            Opt.Threshold,numel(dspec),BasisSize);
        else
          thisspec = ones(size(omega));
          logmsg(0,'  Tridiagonalization did not converge to within %g after %d steps!\n  Increase Options.LLMK (current settings [%d,%d,%d,%d])',...
            Opt.Threshold,BasisSize,Opt.LLMK');
        end
        
      case '\' % MATLAB backslash solver for sparse linear system
        
        I = speye(size(L));
        rho0 = StartVector;
        if explicitFieldSweep
          Q = L - 1i*omega*I;
          thisspec(iB) = rho0'*(Q\rho0);
        else
          for iOmega = 1:numel(omega)
            Q = L - 1i*omega(iOmega)*I;
            thisspec(iOmega) = rho0'*(Q\rho0);
          end
        end
        
      case 'E' % eigenvalue method (sum of Lorentzians)
        % see e.g. G. Binsch, J. Am. Chem. Soc. 91, 1304 (1969)
        L = full(L);
        [U,Lam] = eig(L);
        Lam = diag(Lam);
        rho0 = StartVector;
        Amplitude = (rho0'*U).'.*(U\rho0);
        thisspec = 0;
        for iPeak = 1:numel(Amplitude)
          thisspec = thisspec + Amplitude(iPeak)./(Lam(iPeak)-1i*omega);
        end
        
      case 'B' % bi-conjugate gradients stabilized
        I = speye(size(L));
        rho0 = StartVector;
        for iOmega = 1:numel(omega)
          Q = L - 1i*omega(iOmega)*I;
          [u,~] = bicgstab(Q,rho0,Opt.Threshold,nDim);
          thisspec(iOmega) = rho0'*u;
        end
        
      case 'C' % conjugate gradients
        CGshift = 1e-6 + 1e-6i;
        [~,alpha,beta,err,stepsDone] = chili_conjgrad(L,StartVector,CGshift);        
        logmsg(1,'  step %d/%d: CG converged to within %g',...
          stepsDone,BasisSize,err);
        thisspec = chili_contfracspec(-1i*omega,alpha,beta);
        
    end
    
  end % field loop
  
  spec = spec + thisspec*Weights(iOri);
  
end % orientation loop

spec = spec/scale;

% Rescale to match rigid-limit chili intensities to pepper intensities
spec = spec*1e10;
spec = spec/2; % since chili uses normalized Sx and pepper uses unnormalized Sx
% (works only for S=1/2)

if Opt.highField
  spec = spec/2;
end

if FrequencySweep
  spec = spec*1e3;
end

% Save structure with internal data to workspace for diagnostics
if saveDiagnostics
  assignin('base',Opt.Diagnostics,diagnostics);
end

%===============================================================================



%===============================================================================
% Phasing
%===============================================================================
spec = real(exp(1i*Exp.mwPhase)*spec);


%===============================================================================
% Post-convolution
%===============================================================================
if doPostConvolution
  logmsg(1,'Postconvolution...');
  
  % Spin system with shf nuclei only
  pcidx = Opt.PostConvNucs;
  pcSys.g = mean(fullSys.g(:));
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
  
  EasySpinLogLevel = Opt.Verbosity; % re-set it, since garlic clears it
  
  % Convolute SLE spectrum with isotropic spectrum
  spec = conv(spec,spec_pc,'same');
end


%===============================================================================
% Basis set analysis
%===============================================================================
if ~isfield(Opt,'BasisAnalysis')
  Opt.BasisAnalysis = false;
end
if Opt.BasisAnalysis
  logmsg(1,'-------------------------------------------------------------------');
  logmsg(1,'Basis set analysis');
  omega_ = linspace(omega(1),omega(end),12);
  u_sum = 0;
  for iOmega = 1:numel(omega_)
    u = bicgstab(L+omega_(iOmega)*speye(size(L)),StartVector,1e-7,180);
    u_sum = u_sum + abs(u)/abs(StartVector'*u);
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


%===============================================================================
% Temperature: include Boltzmann equilibrium polarization
%===============================================================================
if isfinite(Exp.Temperature)
  if FieldSweep
     DeltaE = planck*Exp.mwFreq*1e9; % joule
  else
     DeltaE = bmagn*gavg*Exp.Field*1e-3; % joule
  end
  e = exp(-DeltaE/boltzm/Exp.Temperature);
  Population = [1 e];
  Population = Population/sum(Population);
  Polarization = Population(1) - Population(2);
  spec = spec*Polarization;
end


%===============================================================================
% Convolutional broadening
%===============================================================================
% Convolution with Gaussian only. Lorentzian broadening is already
% included in the slow-motion simulation via T2.
fwhmG = Sys.lw(1);
if fwhmG>0 && ConvolutionBroadening
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


%===============================================================================
% Field modulation, or derivatives
%===============================================================================
logmsg(1,'Modulation');

% Pseudo field modulation
if FieldSweep
  if Exp.ModAmp>0
    logmsg(1,'  pseudo field modulation with %g mT, harmonic %d',Exp.ModAmp,Exp.ModHarmonic);
    outspec = fieldmod(xAxis,outspec,Exp.ModAmp,Exp.ModHarmonic);
  else
    logmsg(1,'  no pseudo field modulation');
  end
else
  logmsg(1,'  no modulation');
end

% Derivatives
if Exp.DerivHarmonic>0
  logmsg(1,'  harmonic %d: using differentiation',Exp.DerivHarmonic);
  dx = xAxis(2)-xAxis(1);
  for h = 1:Exp.DerivHarmonic
    dspec = diff(outspec,[],2)/dx;
    outspec = (dspec(:,[1 1:end]) + dspec(:,[1:end end]))/2;
  end
else
  logmsg(1,'  no derivative');
end


%===============================================================================
%  Final processing
%===============================================================================

switch nargout
  case 0
    plotresults(xAxis,outspec,Exp,~FieldSweep);
  case 1
    varargout = {outspec};
  case 2
    varargout = {xAxis,outspec};
end
%===============================================================================


logmsg(1,'-------------------------------------------------------------------');

clear global EasySpinLogLevel

return
%===============================================================================
%===============================================================================
%===============================================================================



%===============================================================================
function Basis = processbasis(Basis,maxPotentialK,I,Symmetry)

nNuclei = numel(I);

nobetatilts = Symmetry.nobetatilts;
tensorsCollinear = Symmetry.tensorsCollinear;
axialSystem = Symmetry.axialSystem;

% Spatial basis parameters: evenLmax oddLmax Kmax Mmax jKmin evenK
%-------------------------------------------------------------------------------
Basis.evenLmax = Basis.LLMK(1);
Basis.oddLmax = Basis.LLMK(2);
Basis.Mmax = Basis.LLMK(3);
Basis.Kmax = Basis.LLMK(4);

% Set jKmin = +1 if tensorial coefficients are all real. This is the
% case when all tensors (g and A) are collinear and tilted relative to
% the diffusion tensor by the angles (0,beta,0).
if isempty(Basis.jKmin)
  Basis.jKmin = -1;
  if tensorsCollinear
    Basis.jKmin = +1;
  end
end

% Use only even K if there is no beta tilt (i.e. all +1 and -1 spherical
% tensor components in the Hamiltonian are zero).
if isempty(Basis.evenK)
  Basis.evenK = nobetatilts;
end

%{
% Use only even L values (oddLmax=0) and no K values (Kmax=0)
% in case of axial magnetic tensors, axial potential, 
% and no magnetic/diffusion tilt
if axialSystem && Basis.evenK && (isempty(maxPotentialK) || (maxPotentialK==0))
  Basis.oddLmax = 0;
  Basis.Kmax = 0;
end
%}

% Spin basis truncation parameters: pSmin, pImax, pImaxall
%-------------------------------------------------------------------------------

% pImax (maximum nuclear coherence order, for each nucleus)
if nNuclei==0
  Basis.pImax = 0;
  Basis.pImaxall = 0;
else
  if isempty(Basis.pImax)
    Basis.pImax = 2*I;
  end
  Basis.pImax = min(Basis.pImax,2*I);
  if isempty(Basis.pImaxall)
    Basis.pImaxall = sum(Basis.pImax);
  end
  Basis.pImaxall = min(Basis.pImaxall,sum(2*I));
end

% Set fields for fast two-nuclei code
if nNuclei==2
  Basis.pI1max = Basis.pImax(1);
  Basis.pI2max = Basis.pImax(2);
end

% Assemble output array of basis set parameters for chili_lm
%-------------------------------------------------------------------------------
deltaK = Basis.evenK+1;
Basis.v = [...
  Basis.evenLmax Basis.oddLmax Basis.Kmax Basis.Mmax, ...
  Basis.jKmin Basis.pSmin deltaK ...
  Basis.MpSymm ...
  Basis.pImax];

return


%===============================================================================
function [Dyn,err] = processdynamics(Dyn,FieldSweep)

err = '';

% Diffusion tensor, correlation time
%-------------------------------------------------------------------------------
% convert everything (tcorr, logcorr, logDiff, Diff) to R
if isfield(Dyn,'Diff')
  % Diff given
  if any(Dyn.Diff<0)
    error('Sys.Diff cannot be negative.');
  end
  Dyn.R = Dyn.Diff;
elseif isfield(Dyn,'logDiff')
  Dyn.R = 10.^Dyn.logDiff;
elseif isfield(Dyn,'tcorr')
  Dyn.R = 1/6./Dyn.tcorr;
elseif isfield(Dyn,'logtcorr')
  if Dyn.logtcorr>=0, error('Sys.logtcorr must be negative.'); end
  Dyn.R = 1/6./10.^Dyn.logtcorr;
else
  err = sprintf('You must specify a rotational correlation time or a diffusion tensor\n(Sys.tcorr, Sys.logtcorr, Sys.Diff or Sys.logDiff).');
  return
end

% check values
if any(Dyn.R<0)
  error('Negative diffusion rate or correlation times are not possible.');
elseif any(Dyn.R>1e12)
  fprintf('Diffusion rate very fast. Simulation might not converge.\n');
elseif any(Dyn.R<1e3)
  fprintf('Diffusion rate very slow. Simulation might not converge.\n');
end

% expand to rhombic tensor
switch numel(Dyn.R)
  case 1, Dyn.R = Dyn.R([1 1 1]);
  case 2, Dyn.R= Dyn.R([1 1 2]);
  case 3 % Diff already rhombic
  otherwise
    err = 'Sys.Diff must have 1, 2 or 3 elements (isotropic, axial, rhombic).';
    return
end

% Linewidth
%-------------------------------------------------------------------------------
if isfield(Dyn,'lw')
  if numel(Dyn.lw)>1
    if FieldSweep
      LorentzFWHM = Dyn.lw(2)*28 * 1e6; % mT -> MHz -> Hz (g = 2.0006)
    else
      LorentzFWHM = Dyn.lw(2)*1e6; % MHz -> Hz
    end
  else
    LorentzFWHM = 0;
  end
  if LorentzFWHM~=0
    % Lorentzian T2 from FWHM in freq domain 1/T2 = pi*FWHM
    Dyn.T2 = 1/LorentzFWHM/pi;
  else
    Dyn.T2 = inf;
  end
end

return


%===============================================================================
function plotresults(xAxis,spec,Exp,FrequencySweep)

cla
if FrequencySweep
  if xAxis(end)<1
    plot(xAxis*1e3,spec);
    xlabel('frequency (MHz)');
  else
    plot(xAxis,spec);
    xlabel('frequency (GHz)');
  end
  title(sprintf('%0.8g mT',Exp.Field));
else
  if xAxis(end)<10000
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

return
