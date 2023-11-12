% pepper  Computation of powder cw EPR spectra 
%
%   pepper(Sys,Exp)
%   pepper(Sys,Exp,Opt)
%   spec = pepper(...)
%   [x,spec] = pepper(...)
%   [x,spec,info] = pepper(...)
%
%   Calculates field-swept and frequency-swept cw EPR spectra.
%
%   Input:
%    Sys: parameters of the paramagnetic system
%      S, g, Nucs, A, Q, D, ee,
%      gFrame, AFrame, QFrame, DFrame, eeFrame
%      lw, lwpp
%      HStrain, gStrain, AStrain, DStrain
%      B2, B4, B6 etc.
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
%      SampleRotation      3-element array of Euler angles (in radians) for sample rotation
%      SampleFrame         3-element array of Euler angles (in radians) for sample/crystal orientations
%      CrystalSymmetry     crystal symmetry (space group etc.)
%      MolFrame            Euler angles (in radians) for molecular frame orientation
%      Mode                excitation mode: 'perpendicular', 'parallel', {k_tilt alpha_pol}
%      Ordering            coefficient for non-isotropic orientational distribution
%    Opt: computational options
%      Method              'matrix', 'perturb1', 'perturb2'='perturb'
%      separate            '', 'components', 'transitions', 'sites', 'orientations'
%      Verbosity           0, 1, 2
%      GridSize            grid size;  N1, [N1 Ninterp]
%      Transitions, Threshold
%      GridSymmetry, GridFrame,
%      Intensity, Freq2Field, Sites
%
%   Output:
%    x        field axis (in mT) or frequency axis (in GHz)
%    spec     spectrum
%    info     structure with details of the calculation
%
%   If no output argument is given, the simulated spectrum is plotted.

function varargout = pepper(Sys,Exp,Opt)

if nargin==0, help(mfilename); return; end

% Check expiry date
error(eschecker);

% Check Matlab version
error(chkmlver);

% Get time for performance report at the end.
StartTime = clock;

% Input argument scanning, get display level and prompt
%=======================================================================
% Guard against wrong number of input or output arguments.
if nargin<1, error('Please supply a spin system as first input argument.'); end
if nargin<2, error('Please supply experimental parameters as second input argument.'); end
if nargin>3, error('Too many input arguments, the maximum is three.'); end
if nargout>3, error('Too many output arguments.'); end

% Initialize options structure to zero if not given.
if nargin<3, Opt = struct; end
if isempty(Opt), Opt = struct; end

if ~isstruct(Sys) && ~iscell(Sys)
  error('The first input (Sys) must be a structure or a cell array of structures.');
end
if ~isstruct(Exp)
  error('The second input (Exp) must be a structure.');
end
if ~isstruct(Opt)
  error('The third input (Opt) must be a structure.');
end

% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;


%==================================================================
% Loop over species and isotopologues
%==================================================================
Exp.FrequencySweep = ~isfield(Exp,'mwFreq') & isfield(Exp,'Field');

if Exp.FrequencySweep
  SweepAutoRange = (~isfield(Exp,'mwRange') || isempty(Exp.mwRange)) && ...
    (~isfield(Exp,'mwCenterSweep') || isempty(Exp.mwCenterSweep));
else
  SweepAutoRange = (~isfield(Exp,'Range') || isempty(Exp.Range)) && ...
    (~isfield(Exp,'CenterSweep') || isempty(Exp.CenterSweep));
end

if ~isfield(Opt,'IsoCutoff'), Opt.IsoCutoff = 1e-3; end


% Process Opt.separate
if ~isfield(Opt,'separate'), Opt.separate = ''; end
[separateOutput,err] = parseoption(Opt,'separate',{'','components','transitions','orientations','sites'});
error(err);
summedOutput = separateOutput==1;
separateComponentSpectra = separateOutput==2;
separateTransitionSpectra = separateOutput==3;
separateOrientationSpectra = separateOutput==4;
separateSiteSpectra = separateOutput==5;
if separateTransitionSpectra || separateSiteSpectra || separateOrientationSpectra
  separateComponentSpectra = true;
end

singleIsotopologue = isfield(Sys,'singleiso') && Sys.singleiso;
if ~singleIsotopologue
  
  thirdOutput = nargout>=3;
  [xAxis,spec,info] = compisoloop(@pepper,Sys,Exp,Opt,SweepAutoRange,thirdOutput,separateComponentSpectra);
  
  % Output and plotting
  switch nargout
    case 0
      cwepr_plot(xAxis,spec,Exp);
    case 1
      varargout = {spec};
    case 2
      varargout = {xAxis,spec};
    case 3
      varargout = {xAxis,spec,info};
  end
  return
end
%==================================================================


%==================================================================
% Single-isotopologue spectrum
%==================================================================

% Now we can start simulating the spectrum
logmsg(1,['=begin=pepper=====' char(datetime) '=================']);
logmsg(2,'  log level %d',EasySpinLogLevel);
logmsg(1,'-general-----------------------------------------------');

%=======================================================================
% Spin system structure
%=======================================================================

[Sys,err] = validatespinsys(Sys);
error(err);

if any(Sys.n>1)
  error('pepper does not support sets of equivalent nuclei (Sys.n>1).');
end

StrainWidths = any([Sys.HStrain(:); Sys.DStrain(:); Sys.gStrain(:); Sys.AStrain(:)]>0);
ConvolutionBroadening = any(Sys.lw>0);

logmsg(1,'  system with %d spin(s) and %d states',numel(spinvec(Sys)),hsdim(Sys));
if StrainWidths, logmsg(1,'  strain widths given'); end
%=======================================================================



%=======================================================================
% Experiment structure, contains experimental settings
%=======================================================================

% Documented fields and their defaults (mandatory parameters are set to NaN)
%DefaultExp.mwFreq = NaN; % for field sweeps
%DefaultExp.Field = NaN; % for frequency sweeps
DefaultExp.CenterSweep = NaN;
DefaultExp.Range = NaN;
DefaultExp.mwCenterSweep = NaN;
DefaultExp.mwRange = NaN;
DefaultExp.nPoints = 1024;
DefaultExp.Temperature = NaN;
DefaultExp.Harmonic = NaN;
DefaultExp.mwMode = 'perpendicular';
DefaultExp.Ordering = [];
DefaultExp.ModAmp = 0;
DefaultExp.mwPhase = 0;
DefaultExp.lightBeam = '';  % no photoexcitation

DefaultExp.SampleRotation = [];
DefaultExp.SampleFrame = [0 0 0];
DefaultExp.CrystalSymmetry = '';
DefaultExp.MolFrame = [];

Exp = adddefaults(Exp,DefaultExp);

% Check for obsolete fields
if isfield(Exp,'Orientations')
  error('Exp.Orientations is no longer supported, use Exp.SampleFrame instead.');
end
if isfield(Exp,'CrystalOrientation')
  error('Exp.CrystalOrientation is no longer supported, use Exp.SampleFrame instead.');
end

% Check microwave frequency and static field
if ~isfield(Exp,'mwFreq') || isempty(Exp.mwFreq)
  if ~isfield(Exp,'Field')
    error('Please supply either the microwave frequency in Exp.mwFreq (for field sweeps) or the magnetic field in Exp.Field (for frequency sweeps).');
  end
  FieldSweep = false;
else
  if isfield(Exp,'Field') && ~isempty(Exp.Field)
    error('Give either Exp.mwFreq (for a field sweep) or Exp.Field (for a frequency sweep), but not both.');
  end
  FieldSweep = true;
end

if ~FieldSweep
  Sys.lw = Sys.lw/1e3;  % MHz -> GHz
end

if FieldSweep
  if numel(Exp.mwFreq)~=1 || any(Exp.mwFreq<=0) || ~isreal(Exp.mwFreq)
    error('Uninterpretable microwave frequency in Exp.mwFreq.');
  end
  logmsg(1,'  field sweep, mw frequency %0.8g GHz',Exp.mwFreq);
else
  if numel(Exp.Field)~=1 || ~isreal(Exp.Field)
    error('Uninterpretable magnetic field in Exp.Field.');
  end
  logmsg(1,'  frequency sweep, magnetic field %0.8g mT',Exp.Field);
end

% Automatic field range determination
if FieldSweep
  if all(isnan(Exp.CenterSweep)) && all(isnan(Exp.Range))
    if numel(Sys.S)==1 && (Sys.S==1/2) && ~any(Sys.L(:))
      logmsg(1,'  automatic determination of sweep range');
      I = nucspin(Sys.Nucs).';
      if ~isempty(I)
        if Sys.fullA
          Amax = max(abs(Sys.A),[],2);
          Amax = max(reshape(Amax,3,[])).';
        else
          Amax = max(abs(Sys.A),[],2);
        end
      else
        Amax = 0;
      end
      hf = sum(I.*Amax)*1e6; % MHz -> Hz
      if Sys.fullg
        for k = 1:Sys.nElectrons
          g(:,k) = eig(Sys.g((1:3)+(k-1)*3,:));  %#ok
        end
        g = g(:);
      else
        g = Sys.g(:);
      end
      gmax = max(g);
      gmin = min(g);
      minB = planck*(Exp.mwFreq*1e9 - hf)/bmagn/gmax/1e-3; % mT
      maxB = planck*(Exp.mwFreq*1e9 + hf)/bmagn/gmin/1e-3; % mT
      Center = (maxB+minB)/2; % mT
      Sweep = maxB-minB+3*max(Sys.lw); % mT
      if Sweep==0, Sweep = 5*max(Sys.lw); end
      if Sweep==0, Sweep = 10; end
      Stretch = 1.25;
      Exp.CenterSweep = [Center, Stretch*Sweep];
    else
      error(sprintf('Cannot automatically determine field range.\nPlease provide either Exp.CenterSweep or Exp.Range.'));
    end
  end
else
  % Automatic range for frequency sweep is done later.
end

% Check both CenterSweep and Range, prefer CenterSweep
if FieldSweep
  if ~isnan(Exp.CenterSweep)
    Exp.Range = Exp.CenterSweep(1) + [-1 1]*Exp.CenterSweep(2)/2;
    if Exp.Range(1)<0
      error('Lower field limit from Exp.CenterSweep cannot be negative.');
    end
  end
  if isfield(Exp,'Range') && all(~isnan(Exp.Range))
    if any(diff(Exp.Range)<=0) || any(~isfinite(Exp.Range)) || ...
        ~isreal(Exp.Range)
      error('Exp.Range is not valid!');
    end
    if any(Exp.Range<0)
      error('Negative magnetic fields in Exp.Range are not possible.');
    end
  end
else
  if ~isnan(Exp.mwCenterSweep)
    Exp.mwRange = Exp.mwCenterSweep(1) + [-1 1]*Exp.mwCenterSweep(2)/2;
    Exp.mwRange = max(Exp.mwRange,0);
  end
  if isfield(Exp,'mwRange') && all(~isnan(Exp.mwRange))
    if diff(Exp.mwRange)<=0 || any(~isfinite(Exp.mwRange)) || ...
        any(~isreal(Exp.mwRange))
      error('Exp.mwRange is not valid!');
    end
    if any(Exp.mwRange<0)
      error('Exp.mwRange cannot be negative.');
    end
  end
end


% Number of points
if any(~isreal(Exp.nPoints)) || numel(Exp.nPoints)>1 || (Exp.nPoints<2)
  error('Problem with Exp.nPoints. Needs to be a number not smaller than 2.')
end

if FieldSweep
  logmsg(1,'  frequency %g GHz, field range [%g %g] mT, %d points',...
    Exp.mwFreq,Exp.Range(1),Exp.Range(2),Exp.nPoints);
else
  if ~SweepAutoRange
    logmsg(1,'  field %g mT, frequency range [%g %g] GHz, %d points',...
      Exp.Field,Exp.mwRange(1),Exp.mwRange(2),Exp.nPoints);
  else
    logmsg(1,'  field %g mT, automatic frequency range, %d points',...
      Exp.Field,Exp.nPoints);
  end
end

% Detection harmonic
autoHarmonic = ~isfield(Exp,'Harmonic') || isempty(Exp.Harmonic) || isnan(Exp.Harmonic);
noBroadening = ~StrainWidths && ~ConvolutionBroadening;
if autoHarmonic
  if FieldSweep && ~noBroadening
    if noBroadening
      Exp.Harmonic = 0;
    else
      Exp.Harmonic = 1;
    end
  else
    Exp.Harmonic = 0;
  end
end
if ~any(Exp.Harmonic==[-1,0,1,2])
  error('Exp.Harmonic must be either 0, 1 or 2.');
end
if noBroadening && (Exp.Harmonic~=0)
  error('\n  No broadening given. Cannot compute spectrum with Exp.Harmonic=%d.\n',Exp.Harmonic);
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

% Microwave phase
if ~FieldSweep
  % flip dispersion lineshape depending on field or freq sweep
  Exp.mwPhase = -Exp.mwPhase;
end


% Resonator mode
if ischar(Exp.mwMode) && ~isempty(Exp.mwMode)
  if strcmp(Exp.mwMode,'perpendicular')
  elseif strcmp(Exp.mwMode,'parallel')
  else
    error('Exp.mwMode must be either ''perpendicular'' or ''parallel''.');
  end
end
logmsg(1,'  harmonic %d, %s mode',Exp.Harmonic,Exp.mwMode);


% Detect sample type (powder/partially ordered vs. crystal)
PowderSimulation = ~isempty(Exp.Ordering) ||(isempty(Exp.MolFrame) && isempty(Exp.CrystalSymmetry));
if PowderSimulation
  if ~isempty(Exp.MolFrame)
    error('Exp.Ordering cannot be used simultaneously with Exp.MolFrame.');
  end
  if ~isempty(Exp.CrystalSymmetry)
    error('Exp.Ordering cannot be used simultaneously with Exp.CrystalSymmetry.');
  end
else
  if isempty(Exp.MolFrame), Exp.MolFrame = [0 0 0]; end
  if isempty(Exp.CrystalSymmetry), Exp.CrystalSymmetry = 'P1'; end
end
Exp.PowderSimulation = PowderSimulation;  % for communication with resf*

% Process Exp.Ordering
if ~isempty(Exp.Ordering)
  if isnumeric(Exp.Ordering) && numel(Exp.Ordering)==1 && isreal(Exp.Ordering)
    lambda = Exp.Ordering;
    Exp.Ordering = @(beta) exp(lambda*plegendre(2,0,cos(beta)));
    logmsg(1,'  partial ordering (built-in function, coefficient = %g)',lambda);
  elseif isa(Exp.Ordering,'function_handle')
    logmsg(1,'  partial ordering (user-supplied function)');
  else
    error('Exp.Ordering must be either a single number or a function handle.');
  end
  if nargin(Exp.Ordering)==1
    Exp.Ordering = @(beta,gamma) Exp.Ordering(beta).*ones(size(gamma));
  elseif nargin(Exp.Ordering)>2
    logmsg(1,'  Ordering function in Exp.Ordering must take 1 argument (beta) or 2 arguments (beta,gamma).');
  end
end

% Temperature and non-equilibrium populations
nonEquiPops = isfield(Sys,'initState') && ~isempty(Sys.initState);
if nonEquiPops
  msg = '  user-specified non-equilibrium state';
else
  if numel(Exp.Temperature)~=1
    error('If given, Exp.Temperature must be a single number.');
  end
  if isfinite(Exp.Temperature)
    msg = sprintf('  temperature %g K',Exp.Temperature);
  else
    msg = '  no temperature';
  end
end
logmsg(1,msg);

[Exp.R_sample,rotateSample] = p_samplerotmatrix(Exp.SampleRotation);

%=======================================================================




%=======================================================================
% Options structure
%=======================================================================
% Documented options, resfields (for the defaults, see resfields).
%DefaultOpt.Transitions = []; % resfields 
%DefaultOpt.Threshold = 1e-3; % resfields
%DefaultOpt.Intensity = 1; % resfields
%DefaultOpt.Freq2Field = 1; % resfields

% Obsolete fields, pepper
obsoleteOptions = {'Convolution','Width'};
for k = 1:numel(obsoleteOptions)
  if isfield(Opt,obsoleteOptions{k})
    error('Options.%s is obsolete. Please remove from code!',obsoleteOptions{k});
  end
end

if isfield(Opt,'nKnots')
  error('Options.nKnots is obsolete. Use Options.GridSize instead, e.g. Options.GridSize = 91.');
end
if isfield(Opt,'Symmetry')
  error('Options.Symmetry is obsolete. Use Options.GridSymmetry instead, e.g. Options.GridSymmetry = ''D2h''.');
end
if isfield(Opt,'nSpline')
  error('Options.nSpline is obsolete. Use a second number in Options.GridSize instead, e.g. Options.GridSize = [19 5] for Options.nSpline = 5.');
end
if isfield(Opt,'Perturb')
  error('Options.Perturb is obsolete. Use Opt.Method=''perturb'' or Opt.Method=''hybrid'' instead.');
end
if isfield(Opt,'Output')
  error('Options.Output is obsolete. Use Options.separate instead.');
end

% Documented fields, pepper
DefaultOpt.Verbosity = 0;
DefaultOpt.GridSymmetry = '';
DefaultOpt.GridFrame = [];
DefaultOpt.separateOutput = '';
DefaultOpt.Method = 'matrix'; % 'matrix', 'eig', 'perturb1', 'perturb2'='perturb' 

% Undocumented fields, pepper
%DefaultOptions.TPSGridSize = 4;      % resfields
%DefaultOptions.TPSGridSymm = 'D2h';  % resfields
%DefaultOpt.Weak = [];                % resfields
DefaultOpt.Smoothing = 2;
DefaultOpt.GridSizeMinimum = 10;
DefaultOpt.Debug = 0;
DefaultOpt.Intensity = 'on';
DefaultOpt.BruteForce = 0;
DefaultOpt.ImmediateBinning = 0;
DefaultOpt.PaddingMultiplier = 3; % for padding before convolution

GridSizeMatrix = [19 4];
GridSizePerturb = [19 4];

Opt = adddefaults(Opt,DefaultOpt);

if FieldSweep
  [Method,err] = parseoption(Opt,'Method',{'eig','matrix','perturb','perturb1','perturb2','hybrid'});
  error(err);
else
  [Method,err] = parseoption(Opt,'Method',{'matrix','perturb','perturb1','perturb2','hybrid'});
  error(err);
  Method = Method + 10;
end

useEigenFields = Method==1;
usePerturbationTheory = any(Method==[3 4 5 12 13 14]);
if Opt.ImmediateBinning && ~usePerturbationTheory
  error('Opt.ImmediateBinning works only with perturbation theory.');
end
if usePerturbationTheory && nonEquiPops
  error('Perturbation theory not available for systems with non-equilibrium populations.');
end

if ~isfield(Opt,'GridSize')
  if usePerturbationTheory
    Opt.GridSize = GridSizePerturb;
  else
    Opt.GridSize = GridSizeMatrix;
  end
end

if Opt.GridSize(1)<Opt.GridSizeMinimum
  error('Options.GridSize must not be less than %d.',Opt.GridSizeMinimum);
end
if numel(Opt.GridSize)<2
  if usePerturbationTheory
    Opt.GridSize(2) = GridSizePerturb(2);
  else
    Opt.GridSize(2) = GridSizeMatrix(2);
  end
end

% Some compatibility checks for separate spectra output (Opt.separate)
if PowderSimulation
  if separateOrientationSpectra
    error(sprintf('\nCannot return separate orientations for powder spectra (Opt.separate=''orientations'').\nUse other setting for Opt.separate.\n'));
  end
  if separateSiteSpectra
    error(sprintf('\nCannot return separate sites for powder spectra (Opt.separate=''sites'').\nUse other setting for Opt.separate.\n'));
  end
else
  if separateTransitionSpectra
    error(sprintf('\n  Cannot return separate transitions for crystal spectra (Opt.separate=''transitions'').\n  Use other setting for Opt.separate.\n'));
  end
end
if Opt.ImmediateBinning && ~summedSpectra
  error('When using Opt.ImmediateBinning, only Opt.separate='''' is possible');
end

% Parse string options
anisotropicIntensities = parseoption(Opt,'Intensity',{'off','on'}) - 1;
Opt.Intensity = anisotropicIntensities;

% Set up grid etc.
[Exp,Opt] = p_symandgrid(Sys,Exp,Opt);
nOrientations = size(Exp.SampleFrame,1);

% Fold orientational distribution function into grid region
if ~isempty(Exp.Ordering)
  orifun = foldoridist(Exp.Ordering,Opt.GridSymmetry);
end

%=======================================================================
%=======================================================================
%                   PEAK DATA COMPUTATION
%=======================================================================
%=======================================================================

logmsg(1,'-resonances--------------------------------------------');
MethodMsg{1} = 'field sweep, eigenfields (Liouville space)';
MethodMsg{2} = 'field sweep, adaptive segmentation (state space)';
MethodMsg{3} = 'field sweep, second-order perturbation theory';
MethodMsg{4} = 'field sweep, first-order perturbation theory';
MethodMsg{5} = 'field sweep, second-order perturbation theory';
MethodMsg{6} = 'field sweep, hybrid (matrix diagonalization for electron spin, perturbation for nuclei)';
MethodMsg{11} = 'frequency sweep, matrix diagonalization';
MethodMsg{12} = 'frequency sweep, second-order perturbation theory';
MethodMsg{13} = 'frequency sweep, first-order perturbation theory';
MethodMsg{14} = 'frequency sweep, second-order perturbation theory';
MethodMsg{15} = 'frequency sweep, hybrid (matrix diagonalization for electron spin, perturbation for nuclei)';
logmsg(1,'  method: %s',MethodMsg{Method});

if FieldSweep
  % Field sweeps
  %---------------------------------------------------------------------------
  if Method==1
    
    % Eigenfield equation
    %------------------------------------------------------------------------
    anisotropicWidths = false;
    if StrainWidths
      logmsg(-inf,'WARNING: Options.Method: eigenfields method -> strains are ignored!');
    end
    
    Exp1 = Exp;
    Exp1.Range = [0 1e8];
    
    logmsg(2,'  -entering resfields_eig----------------------------------');
    [Pdat,Idat] = resfields_eig(Sys,Exp1,Opt);
    logmsg(2,'  -exiting resfields_eig-----------------------------------');
    Wdat = [];
    %Gdat = [];
    Transitions = [];
    
    if nOrientations==1
      Pdat = {Pdat};
      Idat = {Idat};
    end
    nReson = 0;
    for k = 1:nOrientations
      nReson = nReson + numel(Pdat{k});
    end
    logmsg(1,'  %d resonance in total (%g per orientation)',nReson,nReson/nOrientations);
    
  else
    
    % Matrix diagonalization and perturbation methods
    %------------------------------------------------------------------------
    
    if ~PowderSimulation
      %if ~isfield(Opt,'Perturb'), Opt.Perturb = 0; end
    end
    
    % Set search range larger than requested field range
    Exp1 = Exp;
    Exp1.SearchRange = Exp1.Range + 0.2*diff(Exp.Range)*[-1 1];
    Exp1.SearchRange(Exp1.SearchRange<0) = 0;
    
    Exp1.AccumWeights = Exp.OriWeights;
    
    logmsg(2,'  -entering resfields*----------------------------------');
    switch Method
      case {2,6} % matrix diagonalization, hybrid
        [Pdat,Idat,Wdat,Transitions] = resfields(Sys,Exp1,Opt);
      case {3,5} % 2nd-order perturbation theory
        Opt.PerturbOrder = 2;
        [Pdat,Idat,Wdat,Transitions,spec] = resfields_perturb(Sys,Exp1,Opt);
      case 4 % 1st-order perturbation theory
        Opt.PerturbOrder = 1;
        [Pdat,Idat,Wdat,Transitions,spec] = resfields_perturb(Sys,Exp1,Opt);
    end
    logmsg(2,'  -exiting resfields*-----------------------------------');
    
    if isempty(Wdat)
      anisotropicWidths = false;
    else
      anisotropicWidths = max(Wdat(:))>0;
    end
    
  end
  
else
  
  % Frequency sweeps
  %--------------------------------------------------------------------------------
  % Matrix diagonalization and perturbation methods
  %------------------------------------------------------------------------
  
  if ~PowderSimulation
    %if ~isfield(Opt,'Perturb'), Opt.Perturb = 0; end
  end
  
  Exp1 = Exp;  
  Exp1.AccumWeights = Exp.OriWeights;
  
  logmsg(2,'  -entering resfreqs*----------------------------------');
  switch Method
    case {11,15} % matrix diagonalization, hybrid
      [Pdat,Idat,Wdat,Transitions] = resfreqs_matrix(Sys,Exp1,Opt);
    case {12,14} % 2nd-order perturbation theory
      Opt.PerturbOrder = 2;
      [Pdat,Idat,Wdat,Transitions,spec] = resfreqs_perturb(Sys,Exp1,Opt);
    case 13 % 1st-order perturbation theory
      Opt.PerturbOrder = 1;
      [Pdat,Idat,Wdat,Transitions,spec] = resfreqs_perturb(Sys,Exp1,Opt);
  end
  logmsg(2,'  -exiting resfreqs*-----------------------------------');
  Pdat = Pdat/1e3; % MHz -> GHz
  Wdat = Wdat/1e3; % MHz -> GHz
  
  if isempty(Wdat)
    anisotropicWidths = false;
  else
    anisotropicWidths = max(Wdat(:))>0;
  end
  
end

if ~anisotropicIntensities
  logmsg(1,'  neglecting amplitude anisotropy, taking amplitude average');
  if ~useEigenFields, Idat = mean(Idat(:))*ones(size(Idat)); end
end

nTransitions = size(Transitions,1);
if ~PowderSimulation && ~isempty(Exp.CrystalSymmetry)
  nSites = numel(Pdat)/nTransitions/nOrientations;
else
  nSites = 1;
end

% Automatic range detemination for frequency sweeps
if ~FieldSweep && SweepAutoRange
  % everything in GHz
  minFreq = min(Pdat(:));
  maxFreq = max(Pdat(:));
  padding = (maxFreq-minFreq)/5;
  if padding==0, padding = 0.1; end
  if anisotropicWidths
    padding = max(padding,5*max(Wdat(:)));
  end
  padding = max(padding,5*sum(Sys.lw));
  minRange = max(0,minFreq-padding);
  maxRange = maxFreq + padding;
  Exp.mwRange = [minRange maxRange]; % GHz
  logmsg(1,'  automatic frequency range [%g %g] GHz',...
    Exp.mwRange(1),Exp.mwRange(2));
end

%=======================================================================
%=======================================================================


% Analysis of the resonance data
%------------------------------------------------------------------------
% If looping fields or out-of-range fields are encountered,
% skip interpolation and projection.
if ~useEigenFields
  NaN_in_Pdat = any(isnan(Pdat),2);
else
  NaN_in_Pdat = 0;
end

if FieldSweep
  if Method~=6
    LoopingTransitionsPresent = size(unique(Transitions,'rows'),1)<size(Transitions,1);
    if LoopingTransitionsPresent && PowderSimulation
      logmsg(0,'** Looping transitions found. Artifacts at coalescence points possible.');
    end
  else
    % hybrid method: Transitions contains replicas of core sys transitions
    LoopingTransitionsPresent = 0;
  end
else
  LoopingTransitionsPresent = 0;
end

if FieldSweep
  SweepRange = Exp.Range;
else
  SweepRange = Exp.mwRange;
end
xAxis = linspace(SweepRange(1),SweepRange(2),Exp.nPoints);
Exp.deltaX = xAxis(2)-xAxis(1);

%=======================================================================
%=======================================================================
%    SPECTRUM CONSTRUCTION (incl. INTERPOLATION and PROJECTION)
%=======================================================================
%=======================================================================
% The position/amplitude/width data from above are
%     (1) interpolated to get more of them, and then
%     (2) projected to obtain the spectrum
% The interpolation and projection algorithms depend
% on the symmetry of the grid.

logmsg(1,'-absorption spectrum construction----------------------');

BruteForceSum = useEigenFields | Opt.BruteForce;
axialGrid = Opt.nOctants==0;
usingGrid = Opt.nOctants>=0;

if Opt.ImmediateBinning

  % Lines have already been binned into the spectral vector on the fly during calculation
  
elseif ~BruteForceSum
  
  % Preparations for interpolation
  %-----------------------------------------------------------------------
  doInterpolation = usingGrid && Opt.GridSize(2)>1;
  if doInterpolation
    % Set an option for the sparse tridiagonal matrix \ solver in global cubic
    % spline interpolation. This function needs some time, so it was taken
    % out of Matlab's original spline() function, which is called many times.
    spparms('autommd',0);
    % Interpolation parameters. 1st char: g global, l linear. 2nd char: order.
    if axialGrid % axial symmetry: 1D interpolation
      if any(NaN_in_Pdat)
        InterpMode = {'L3','L3','L3'};
      else
        InterpMode = {'G3','L3','L3'};
      end
    else % 2D interpolation (no L3 method available)
      if any(NaN_in_Pdat)
        InterpMode = {'L1','L1','L1'};
      else
        InterpMode = {'G3','L1','L1'};
      end
    end
    msg = sprintf('  interpolation (factor %d, method %s/%s/%s)',Opt.GridSize(2),InterpMode{1},InterpMode{2},InterpMode{3});
  else
    Opt.GridSize(2) = 1;
    msg = '  interpolation off';
  end
  nfKnots = (Opt.GridSize(1)-1)*Opt.GridSize(2) + 1;
  logmsg(1,msg);
  
  % Preparations for summation/projection
  %-----------------------------------------------------------------------
  doProjection = usingGrid && ~anisotropicWidths;
  if ~doProjection
    msg = 'summation';
    % Construct Gaussian template lineshape
    x0T = 5e4;  % center
    wT = x0T/2.5;  % width; results in <2e-9 at borders
    xT = 0:2*x0T-1;  % needs to be zero-based and with increment 1 (for lisum1i)
    Template = gaussian(xT,x0T,wT,-1);
    if ~anisotropicWidths
      if Sys.lw(1)>0
        % If present, use convolutional Gaussian for spectrum construction
        thisWid = Sys.lw(1);
        Sys.lw(1) = 0;
        ConvolutionBroadening = Sys.lw(2)>0;
        % In the absence of additional convolutional broadening, use derivative
        % to calculate harmonic
        if ~ConvolutionBroadening
          Exp.DerivHarmonic = Exp.DerivHarmonic+Exp.ConvHarmonic;
          Exp.ConvHarmonic = 0;
        end
      else
        % Summation & Lorentzian only: use Lorentzian template
        x0T = 1e5;
        wT = x0T/20; % 0.0025 at borders for Harmonic = -1
        xT = 0:2*x0T-1;
        Template = lorentzian(xT,x0T,wT,Exp.ConvHarmonic-1,Exp.mwPhase);
        thisWid = Sys.lw(2);
        Sys.lw(2) = 0;
        ConvolutionBroadening = false;
      end
      if ~PowderSimulation
        thisWid = repmat(thisWid,nTransitions,1);
      end
    else
      % thisWid is assigned inside the transition/orientation/site loop
    end
  else
    msg = 'triangle/segment projection';
  end
  logmsg(1,'  %s',msg);

  % Pre-allocation of spectral array
  %-----------------------------------------------------------------------
  if separateOrientationSpectra
    nSpectra = nOrientations;
    msg = '(separate orientations)';
  elseif separateSiteSpectra
    nSpectra = nSites;
    msg = '(separate sites)';
  elseif separateTransitionSpectra
    nSpectra = nTransitions;
    msg = '(separate transitions)';
  else
    nSpectra = 1;
  end
  spec = zeros(nSpectra,Exp.nPoints);
  logmsg(1,'  spectrum array size: %dx%d %s',size(spec,1),size(spec,2),msg);  
  
  % Spectrum construction
  %-----------------------------------------------------------------------
  if ~PowderSimulation
    %=======================================================================
    % Single-crystal spectra
    %=======================================================================
    
    if nTransitions>0
      if ~anisotropicIntensities, thisInt = ones(nTransitions,1); end
      %if ~anisotropicWidths, thisWid = zeros(nTransitions,1); end
      
      iOriSite = 1;  % index into Pdat/Idat/Wdat
      spcidx = 0;  % index into spectral output array
      for iOri = 1:nOrientations
        if separateSiteSpectra, spcidx = 0;
        elseif separateOrientationSpectra, spcidx = spcidx+1;
        else, spcidx = 1;
        end
        for iSite = 1:nSites
          if separateSiteSpectra, spcidx = spcidx + 1; end
          %logmsg(3,'  orientation %d of %d, site %d of %d',iOri,nOrientations,iSite,nSites);

          thisPos = Pdat(:,iOriSite);
          if anisotropicIntensities, thisInt = Idat(:,iOriSite); end
          if anisotropicWidths, thisWid = Wdat(:,iOriSite); end
          
          thisspec = lisum1i(Template,x0T,wT,thisPos,thisInt,thisWid,xAxis);
          thisspec = thisspec/nSites/nOrientations;
          thisspec = (2*pi)*thisspec; % for consistency with powder spectra (factor from integral over chi)
          thisspec = Exp.OriWeights(iOri)*thisspec; % integral over (phi,theta)
          
          if ~separateSiteSpectra && ~separateOrientationSpectra
            spec = spec + thisspec;
          else
            spec(spcidx,:) = spec(spcidx,:) + thisspec;
          end
          
          iOriSite = iOriSite + 1;
        end
      end
    end
    
  elseif ~usingGrid
    
    %=======================================================================
    % Isotropic powder spectra
    %=======================================================================
    
    if ~anisotropicIntensities, thisInt = 1; end
    %if ~anisotropicWidths, thisWid = 0; end
    
    spcidx = 0;
    for iTrans = 1:nTransitions
      %logmsg(3,'  transition %d of %d',iTrans,nTransitions);

      thisPos = Pdat(iTrans,:);
      if anisotropicIntensities, thisInt = Idat(iTrans,:); end
      if anisotropicWidths, thisWid = Wdat(iTrans,:); end
      
      thisspec = lisum1i(Template,x0T,wT,thisPos,thisInt,thisWid,xAxis);
      thisspec = (2*pi)*thisspec; % integral over chi (0..2*pi)
      thisspec = Exp.OriWeights*thisspec; % integral over (phi,theta)
      
      if ~separateTransitionSpectra
        spec = spec + thisspec;
      else
        spcidx = spcidx + 1;
        spec(spcidx,:) = thisspec;
      end
      
    end
    
  else
    
    %=======================================================================
    % Anisotropic powder spectra: interpolation and accumulation/projection
    %=======================================================================
    if axialGrid
      if doInterpolation
        % set up fine interpolation grid
        grid = sphgrid(0,nfKnots);
        fphi = grid.phi;
        fthe = grid.theta;
      else
        fthe = Exp.theta;
      end
      fSegWeights = -diff(cos(fthe))*4*pi; % sum is 4*pi
      
      % Obtain user-supplied orientational distribution weights
      if ~isempty(Exp.Ordering)
        centerTheta = (fthe(1:end-1)+fthe(2:end))/2;
        centerPhi = zeros(1,numel(centerTheta));
        if rotateSample
          v = ang2vec(centerPhi,centerTheta);
          [centerPhi,centerTheta] = vec2ang(Exp.R_sample*v);
        end
        OrderingWeights = orifun(-centerTheta,-centerPhi);
        if any(OrderingWeights<0), error('User-supplied orientation distribution gives negative values.'); end
        if all(OrderingWeights==0), error('User-supplied orientation distribution is all-zero.'); end
        fSegWeights = fSegWeights(:).*OrderingWeights(:);
        fSegWeights = 4*pi/sum(fSegWeights)*fSegWeights;
      end

      logmsg(1,'  total %d segments, %d transitions',numel(fthe)-1,nTransitions);
      
    else % nonaxial grid symmetry
      if doInterpolation
        % set up fine interpolation grid
        [grid,tri] = sphgrid(Opt.GridSymmetry,nfKnots);
        fphi = grid.phi;
        fthe = grid.theta;
      else
        tri = Exp.tri;
        fthe = Exp.theta;
        fphi = Exp.phi;
      end
      idxTri = tri.idx.';
      Areas = tri.areas;

      % Obtain user-supplied orientational distribution weights
      if ~isempty(Exp.Ordering)
        centerTheta = mean(fthe(idxTri));
        centerPhi = mean(fphi(idxTri));
        if rotateSample
          v = ang2vec(centerPhi,centerTheta);
          [centerPhi,centerTheta] = vec2ang(Exp.R_sample*v);
        end
        OrderingWeights = orifun(-centerTheta,-centerPhi);
        if any(OrderingWeights<0), error('User-supplied orientation distribution gives negative values!'); end
        if all(OrderingWeights==0), error('User-supplied orientation distribution is all-zero.'); end
        Areas = Areas(:).*OrderingWeights(:);
        Areas = 4*pi/sum(Areas)*Areas;
      end

      logmsg(1,'  total %d triangles (%d orientations), %d transitions',size(idxTri,2),numel(fthe),nTransitions);
    end
    
    if ~anisotropicIntensities, fInt = ones(size(fthe)); end
    if ~anisotropicWidths, fWid = zeros(size(fthe)); end
    
    minBroadening = inf;
    nBroadenings = 0;
    sumBroadenings = 0;
    spcidx = 0;
    
    for iTrans = 1:nTransitions
      
      % Interpolation
      %------------------------------------------------------
      %LoopTransition = any(isnan(Pdat(iTrans,:)));
      LoopTransition = false;
      interpolateThis = doInterpolation && ~LoopTransition;
      if interpolateThis
        fPos = gridinterp(Pdat(iTrans,:),Opt.GridParams,fphi,fthe,InterpMode{1});
        if anisotropicIntensities
          fInt = gridinterp(Idat(iTrans,:),Opt.GridParams,fphi,fthe,InterpMode{2});
        end
        if anisotropicWidths
          fWid = gridinterp(Wdat(iTrans,:),Opt.GridParams,fphi,fthe,InterpMode{3});
        end
      else
        fPos = Pdat(iTrans,:);
        if anisotropicIntensities, fInt = Idat(iTrans,:); end
        if anisotropicWidths, fWid = Wdat(iTrans,:); end
      end
      
      msg1 = '';
      if ~nonEquiPops && any(fInt(:)<0), msg1 = 'intensities'; end
      if any(fWid(:)<0), msg1 = 'widths'; end
      if ~isempty(msg1)
        error('Negative %s encountered! Please report!',msg1);
      end
      
      % Summation or projection
      %------------------------------------------------------
      projectThis = doProjection && ~LoopTransition;
      if projectThis
        if axialGrid
          thisspec = projectzones(fPos,fInt,fSegWeights,xAxis);
        else
          thisspec = projecttriangles(idxTri,Areas,fPos,fInt,xAxis);
        end
        % minBroadening = ?
      else % do summation
        if axialGrid
          fPosC = (fPos(1:end-1) + fPos(2:end))/2;
          fIntC = fSegWeights(:).'.*(fInt(1:end-1) + fInt(2:end))/2;
          fSpread = abs(fPos(1:end-1) - fPos(2:end));
          fWidM  = (fWid(1:end-1) + fWid(2:end))/2;
          c1 = 1.57246; c2 = 18.6348;
        else
          fPosSorted = sort(fPos(idxTri),1);
          fPosC = mean(fPosSorted,1);
          fIntC = Areas(:).'.*mean(fInt(idxTri),1);
          fSpread = fPosSorted(3,:) - fPosSorted(1,:);
          fWidM = mean(fWid(idxTri),1);
          c1 = 2.8269; c2 = 42.6843;
        end
        Lambda = fWidM./fSpread;
        gam = 1./sqrt(c1*Lambda.^2 + c2*Lambda.^4);
        gam(isinf(gam)) = 0;
        fWidC = fWidM.*(1 + Opt.Smoothing*gam);
        
        thisspec = lisum1i(Template,x0T,wT,fPosC,fIntC,fWidC,xAxis);
        
        minBroadening = min(minBroadening,min(Lambda));
        sumBroadenings = sumBroadenings + sum(Lambda);
        nBroadenings = nBroadenings + numel(Lambda);
      end
      
      thisspec = thisspec*(2*pi); % integral over chi (0..2*pi)
      
      % Accumulate subspectra
      %----------------------------------------------------------
      if ~separateTransitionSpectra
        spec = spec + thisspec;
      else
        spcidx = spcidx + 1;
        spec(spcidx,:) = thisspec;
      end
      
    end % for iTrans
    
    if ~doProjection
      logmsg(1,'  Smoothness: overall %0.4g, worst %0.4g\n   (<0.5: probably bad, 0.5-3: ok, >3: overdone)',sumBroadenings/nBroadenings,minBroadening);
    end
    
  end
  %=======================================================================
  
else % if Opt.ImmediateBinning elseif ~BruteForceSum ...
  
  logmsg(1,'  no interpolation',nOrientations);
  logmsg(1,'  constructing stick spectrum');
  logmsg(1,'  summation over %d orientations',nOrientations);
  spec = zeros(1,Exp.nPoints);
  prefactor = (Exp.nPoints-1)/(SweepRange(2)-SweepRange(1));
  for iOri = 1:nOrientations
    if iscell(Pdat)
      thisP = Pdat{iOri};
      Amplitudes = Idat{iOri};
    else
      thisP = Pdat(:,iOri);
      Amplitudes = Idat(:,iOri);
    end
    idxPos = round(1+prefactor*(thisP-SweepRange(1)));
    outOfRange = (idxPos<1) | (idxPos>Exp.nPoints);
    idxPos(outOfRange) = [];
    Amplitudes(outOfRange) = [];
    if anisotropicIntensities
      spec = spec + full(sparse(1,idxPos,Exp.OriWeights(iOri)*Amplitudes,1,Exp.nPoints));
    else
      spec = spec + full(sparse(1,idxPos,Exp.OriWeights(iOri),1,Exp.nPoints));
    end
  end
  
  spec = spec * (2*pi); % integral over chi (0..2*pi)
  
  spec = spec/Exp.deltaX;
  
end
%=======================================================================





%=======================================================================
%                         Final activities
%=======================================================================
logmsg(1,'-final-------------------------------------------------');

% Combine branches of looping transitions if separate output
%-----------------------------------------------------------------------
if FieldSweep && PowderSimulation
  if ~summedOutput && LoopingTransitionsPresent
    [Transitions,~,idx] = unique(Transitions,'rows');
    nTransitions = size(Transitions,1);
    newspec = zeros(nTransitions,Exp.nPoints);
    for k = 1:length(idx)
      newspec(idx(k),:) = newspec(idx(k),:) + spec(k,:);
    end
    spec = newspec;
    clear newspec;
  end
end


% Convolution with line shape
%-----------------------------------------------------------------------
if ConvolutionBroadening
  logmsg(1,'  harmonic %d: using convolution',Exp.ConvHarmonic);
  fwhmG = Sys.lw(1);
  fwhmL = Sys.lw(2);
  if fwhmL>0
    HarmonicL = Exp.ConvHarmonic;
    mwPhaseL = Exp.mwPhase;
    HarmonicG = 0;
    mwPhaseG = 0;
  else
    HarmonicL = 0;
    mwPhaseL = 0;
    HarmonicG = Exp.ConvHarmonic;
    mwPhaseG = Exp.mwPhase;
  end
  if FieldSweep
    unitstr = 'mT';
  else
    unitstr = 'GHz';
  end
  
  % Add padding to left and right of spectral range to reduce convolution artifacts
  RangePadding = true;
  if RangePadding
    exceedsLowerLimit = any(spec(:,1)~=0);
    exceedsHigherLimit = any(spec(:,end)~=0);
    if exceedsLowerLimit
      if exceedsHigherLimit
        logmsg(0,'** Spectrum exceeds sweep range. Artifacts at lower and upper limits possible.');
      else
        logmsg(0,'** Spectrum exceeds sweep range. Artifacts at lower limit possible.');
      end
    else
      if exceedsHigherLimit
        logmsg(0,'** Spectrum exceeds sweep range. Artifacts at upper limit possible.');
      end
    end
    if  exceedsLowerLimit || exceedsHigherLimit
      nPad = round(max([fwhmG fwhmL])/Exp.deltaX*Opt.PaddingMultiplier);
      spec = [repmat(spec(:,1),1,nPad) spec]; % left padding
      spec = [spec repmat(spec(:,end),1,nPad)]; % right padding
    else
      nPad = 0;
    end
  end

  % Convolution with Lorentzian
  if fwhmL~=0
    if fwhmL>2*Exp.deltaX
      logmsg(1,'  convoluting with Lorentzian, FWHM %g %s, derivative %d',fwhmL,unitstr,HarmonicL);
      if size(spec,1)>1, fwhm = [0 fwhmL]; else, fwhm = fwhmL; end
      spec = convspec(spec,Exp.deltaX,fwhm,HarmonicL,0,mwPhaseL);
    else
      if HarmonicL==0
        % Skip convolution, since it has no effect with such a narrow delta-like Lorentzian.
      else
        error('Lorentzian linewidth (FWHM %g %s) is smaller than 2 increments (2x%g = %g %s) - cannot do convolution.\nIncrease linewidth, or increment number of points in Exp.nPoints.',fwhmL,unitstr,Exp.deltaX,2*Exp.deltaX,unitstr);
      end
    end
  end
  
  % Convolution with Gaussian
  if fwhmG~=0
    if fwhmG>2*Exp.deltaX
      logmsg(1,'  convoluting with Gaussian, FWHM %g %s, derivative %d',fwhmG,unitstr,HarmonicG);
      if size(spec,1)>1, fwhm = [0 fwhmG]; else, fwhm = fwhmG; end
      spec = convspec(spec,Exp.deltaX,fwhm,HarmonicG,1,mwPhaseG);
    else
      if HarmonicG==0
        % Skip convolution, since it has no effect with such a narrow delta-like Gaussian.
      else
        error('Gaussian linewidth (FWHM %g %s) is smaller than 2 increments (2x%g = %g %s) - cannot do convolution.\nIncrease linewidth, increment number of points in Exp.nPoints.',fwhmG,unitstr,Exp.deltaX,2*Exp.deltaX,unitstr);
      end
    end
  end

  % Remove padding
  if RangePadding
    if nPad>0
      spec(:,1:nPad) = [];
      spec(:,Exp.nPoints+1:end) = [];
    end
  end

else
  
  if Exp.DerivHarmonic>0
    logmsg(1,'  harmonic %d: using differentiation',Exp.DerivHarmonic);
    for h = 1:Exp.DerivHarmonic
      dspec = diff(spec,[],2)/Exp.deltaX;
      spec = (dspec(:,[1 1:end]) + dspec(:,[1:end end]))/2;
    end
  else
    logmsg(1,'  harmonic 0: absorption spectrum');
  end

end


% Field modulation
%-----------------------------------------------------------------------
if FieldSweep
  if Exp.ModAmp>0
    logmsg(1,'  applying field modulation');
    if summedOutput
      spec = fieldmod(xAxis,spec,Exp.ModAmp,Exp.ModHarmonic);
    else
      for iSpec = 1:size(spec,1)
        spec(iSpec,:) = fieldmod(xAxis,spec(iSpec,:),Exp.ModAmp,Exp.ModHarmonic);
      end
    end
  else
    % derivatives already included in convolutions etc.
  end
else
  % frequency sweeps: not available
end

% Assign output
%-----------------------------------------------------------------------
switch nargout
  case 1
    varargout = {spec};
  case 2
    varargout = {xAxis,spec};
  case 3
    info.Transitions = Transitions;
    info.nSites = nSites;
    info.nOrientations = nOrientations;
    varargout = {xAxis,spec,info};
end

% Report performance
%-----------------------------------------------------------------------
hmsString = elapsedtime(StartTime,clock);
logmsg(1,['pepper took ' hmsString]);

logmsg(1,'=end=pepper=======%s=================\n',char(datetime));

clear global EasySpinLogLevel

end
