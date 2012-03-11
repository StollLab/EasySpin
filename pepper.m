% pepper  Computation of powder cw EPR spectra 
%
%   pepper(Sys,Exp)
%   pepper(Sys,Exp,Opt)
%   spec = pepper(...)
%   [B,spec] = pepper(...)
%   [B,spec,Trans] = pepper(...)
%
%   Calculates a field-swept cw EPR spectrum.
%
%   Input:
%   - Sys: parameters of the paramagnetic system
%       S, g, Nucs, A, Q, D, ee,
%       gpa, Apa, Qpa, Dpa, eepa
%       lw, lwpp
%       HStrain, gStrain, AStrain, DStrain
%       B20 etc.
%   - Exp: experimental parameters
%       mwFreq        spectrometer frequency, in GHz
%       CenterSweep   magnetic field range [center,sweep], in mT
%       Range         magnetic field range [low,high], in mT
%       nPoints       number of points
%       Temperature   temperature of the sample, by default off (NaN)
%       Harmonic      detection harmonic: 0, 1 (default), 2
%       Mode          resonator mode: 'parallel', 'perpendicular' (default)
%       Orientations  orientations for single-crystal simulations
%       Ordering      coefficient for non-isotropic orientational distribution
%   - Opt: computational options
%       Method        'matrix', 'perturb1', 'perturb2'='perturb'
%       Verbosity, Output,
%       Symmetry, SymmFrame,
%       nKnots, Intensity
%       Transitions, nTransitions, Threshold
%
%   Output:
%   - B:      the field axis, in mT
%   - spec:   the spectrum
%   - Trans:  transitions included in the calculation
%
%   If no output argument is given, the simulated spectrum is plotted.

function varargout = pepper(Sys,Exp,Opt)

if (nargin==0), help(mfilename); return; end

% Get time for performance report at the end.
StartTime = clock;

% Input argument scanning, get display level and prompt
%=======================================================================
% Check Matlab version
VersionErrorStr = chkmlver;
error(VersionErrorStr);

% --------License ------------------------------------------------
LicErr = 'Could not determine license.';
Link = 'epr@eth'; eschecker; error(LicErr); Link; clear Link LicErr
% --------License ------------------------------------------------

% Guard against wrong number of input or output arguments.
if (nargin<1), error('Please supply a spin system as first parameter.'); end
if (nargin<2), error('Please supply experimental parameters as second input argument.'); end
if (nargin>3), error('Too many input arguments, the maximum is three.'); end

if (nargout>4), error('Too many output arguments.'); end

% Initialize options structure to zero if not given.
if (nargin<3), Opt = struct('unused',NaN); end
if isempty(Opt), Opt = struct('unused',NaN); end

if ~isstruct(Sys) && ~iscell(Sys)
  error('Sys must be a structure or a list of structures!');
end
if ~isstruct(Exp)
  error('Exp must be a structure!');
end
if ~isstruct(Opt)
  error('Opt must be a structure!');
end

% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

%==================================================================
% Multiple components

if ~isfield(Sys,'singlecomponent')
  
  if isstruct(Sys), Sys = {Sys}; end
  
  % Error if no field range given for multi-component simulation
  if numel(Sys)>1
    if ~isfield(Exp,'Range') || isempty(Exp.Range)
      if ~isfield(Exp,'CenterSweep') || isempty(Exp.CenterSweep)
        error('For multiple components, please provide a field range in Exp.Range or Exp.CenterSweep.');
      end
    end
  end
  
  % Loop over components
  spec = 0;
  for k = 1:numel(Sys)
    Sys{k}.singlecomponent = 1;
    [xAxis,spec_,transitions] = pepper(Sys{k},Exp,Opt);
    if isfield(Sys{k},'weight')
      spec = spec + spec_*Sys{k}.weight;
    else
      spec = spec + spec_;
    end
  end

  % Output and plotting
  switch nargout
    case 0
      cla
      if (xAxis(2)<10000)
        plot(xAxis,spec);
        xlabel('magnetic field (mT)');
      else
        plot(xAxis/1e3,spec);
        xlabel('magnetic field (T)');
      end
      axis tight
      ylabel('intensity (arb.u.)');
      title(sprintf('%0.8g GHz',Exp.mwFreq));
    case 1, varargout = {spec};
    case 2, varargout = {xAxis,spec};
    case 3, varargout = {xAxis,spec,transitions};
  end
  
  return

end

%==================================================================
% Loop over isotopologues, if necessary. 
if ~isfield(Sys,'singleiso')
  if ~isfield(Opt,'IsoCutoff'), Opt.IsoCutoff = 0.001; end

  Sys.singleiso = 1;
  if ~isfield(Sys,'Nucs'), Sys.Nucs = ''; end
  if ~isfield(Sys,'Abund'), Sys.Abund = []; end
  out = isotopologues(Sys.Nucs,Sys.Abund,[],Opt.IsoCutoff);
  
  if (out.nIso>1)
    AutoRangeRequested = ~isfield(Exp,'Range') & ~isfield(Exp,'CenterSweep');
    if AutoRangeRequested
      error(sprintf('Cannot automatically determine field range.\nPlease manually specify the magnetic field range using Exp.Range or Exp.CenterSweep.'));
    end
  end

  % Loop over isotope combinations
  spec = 0;
  for iIso = 1:out.nIso
    q = out.Nucs{iIso};
    if ~isempty(q)
      Sys.Nucs = q;
      Sys.Ascale = out.Ascale{iIso};
      Sys.Qscale = out.Qscale{iIso};
    else
      Sys.Nucs = [];
    end
    [xAxis,spec_,transitions] = pepper(Sys,Exp,Opt);
    if iIso>1
      if size(spec_,1)~=size(spec,1)
        error('Cannot combine spectra. Set Opt.Output to ''summmed''.');
      end
    end
    spec = spec + out.Abund(iIso)*spec_;
  end

  % Output and plotting
  switch nargout
    case 0
      cla
      if (xAxis(2)<10000)
        plot(xAxis,spec);
        xlabel('magnetic field (mT)');
      else
        plot(xAxis/1e3,spec);
        xlabel('magnetic field (T)');
      end
      axis tight
      ylabel('intensity (arb.u.)');
      title(sprintf('%0.8g GHz',Exp.mwFreq));
    case 1, varargout = {spec};
    case 2, varargout = {xAxis,spec};
    case 3, varargout = {xAxis,spec,transitions};
  end
  return

end
%==================================================================


% Single-isotopologue spectrum
%==================================================================

% Now we can start simulating the spectrum
logmsg(1,['=begin=pepper=====' datestr(now) '=================']);
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
ConvWidth = any(Sys.lw>0);

logmsg(1,'  system with %d spin(s) and %d states',numel(spinvec(Sys)),hsdim(Sys));
if StrainWidths, logmsg(1,'  widths for Gaussian strains given'); end
%=======================================================================



%=======================================================================
% Experiment structure, contains experimental settings
%=======================================================================


% Documented fields and their defaults (mandatory parameters are set to NaN)
DefaultExp.mwFreq = NaN;
DefaultExp.CenterSweep = NaN;
DefaultExp.Range = NaN;
DefaultExp.nPoints = 1024;
DefaultExp.Temperature = NaN;
DefaultExp.Harmonic = 1;
DefaultExp.Mode = 'perpendicular';
DefaultExp.Orientations = [];
DefaultExp.Ordering = [];
DefaultExp.CrystalSymmetry = [];
DefaultExp.ModAmp = 0;
DefaultExp.mwPhase = 0;

Exp = adddefaults(Exp,DefaultExp);

% Check microwave frequency
if ~isfield(Exp,'mwFreq') || any(isnan(Exp.mwFreq))
  error('Please supply the microwave frequency in Exp.mwFreq.');
end
if (numel(Exp.mwFreq)~=1) || any(Exp.mwFreq<=0) || ~isreal(Exp.mwFreq)
  error('Unintelligible microwave frequency in Exp.mwFreq.');
end

% Automatic field range determination
if isnan(Exp.CenterSweep) & isnan(Exp.Range)
  if (Sys.S==1/2)
    logmsg(1,'  automatic determination of magnetic field range');
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
        g(:,k) = eig(Sys.g((1:3)+(k-1)*3,:));
      end
      g = g(:);
    else
      g = Sys.g(:);
    end
    gmax = max(g);
    gmin = min(g);
    minB = planck*(Exp.mwFreq*1e9 - hf)/bmagn/gmax/1e-3;
    maxB = planck*(Exp.mwFreq*1e9 + hf)/bmagn/gmin/1e-3;
    Sweep = maxB-minB;
    Center = (maxB+minB)/2;
    if Sweep==0, Sweep = 5*max(Sys.lw); end
    if Sweep==0, Sweep = 10; end
    Exp.CenterSweep = [Center, 1.25*Sweep];
  else
    error(sprintf('Cannot automatically determine field range.\nPlease given either Exp.CenterSweep or Exp.Range.'));
  end
end

% Check both CenterSweep and Range, prefer CenterSweep
if ~isnan(Exp.CenterSweep)
  Exp.Range = Exp.CenterSweep(1) + [-1 1]*Exp.CenterSweep(2)/2;
  Exp.Range = max(Exp.Range,0);
end

if (diff(Exp.Range)<=0) | ~isfinite(Exp.Range) | ~isreal(Exp.Range) | any(Exp.Range<0)
  %Exp.Range
  %Sys
  error('Exp.Range is not valid!');
end

% Number of points
if any(~isreal(Exp.nPoints)) || numel(Exp.nPoints)>1 || (Exp.nPoints<2)
  error('Problem with Exp.nPoints. Needs to be a number not smaller than 2.')
end

logmsg(1,'  mw %g GHz, range [%g %g] mT, %d points',...
  Exp.mwFreq,Exp.Range(1),Exp.Range(2),Exp.nPoints);


% Detection harmonic
if ~any(Exp.Harmonic==[-1,0,1,2])
  error('Exp.Harmonic must be 0, 1 or 2.');
end
noBroadening = (~StrainWidths) && (~ConvWidth);
if (Exp.Harmonic>0) && noBroadening
  error(['No or zero linewidth/broadening given. Cannot compute spectrum with Exp.Harmonic=%d.\n'...
    'Please specify a line broadening (lwpp, lw, gStrain, AStrain, DStrain).'],Exp.Harmonic);
end

% Resonator mode
if isfield(Exp,'Detection')
  error('Exp.Detection is obsolete. Use Exp.Mode instead.');
end
if strcmp(Exp.Mode,'perpendicular')
  ParallelMode = 0;
elseif strcmp(Exp.Mode,'parallel')
  ParallelMode = 1;
else
  error('Exp.Mode must be either ''perpendicular'' or ''parallel''.');
end
logmsg(1,'  harmonic %d, %s mode',Exp.Harmonic,Exp.Mode);


% Modulation amplitude
if any(Exp.ModAmp<0) || any(isnan(Exp.ModAmp)) || numel(Exp.ModAmp)~=1
  error('Exp.ModAmp must be either a single positive number or zero.');
end
if (Exp.ModAmp>0)
  logmsg(1,'  field modulation, amplitude %g mT',Exp.ModAmp);
  if (Exp.Harmonic<1)
    error('With field modulation (Exp.ModAmp), Exp.Harmonic=0 does not work.');
  end
  Exp.Harmonic = Exp.Harmonic - 1;
end


if isfield(Exp,'Orientation')
  disp('Exp.Orientation given, did you mean Exp.Orientations?');
end
PowderSimulation = isempty(Exp.Orientations);
if ~PowderSimulation
  % Make sure Exp.Orientations is ok
  [n1,n2] = size(Exp.Orientations);
  if ((n2==2)||(n2==3)) && (n1~=2) && (n1~=3)
    Exp.Orientations = Exp.Orientations.';
  end
  [nAngles,nOrientations] = size(Exp.Orientations);
  if (nAngles<2) || (nAngles>3)
    error('Exp.Orientations array has %d rows instead of 2 or 3.',nAngles);
  end
  % Make sure Exp.CrystalSymmetry is ok
  if isempty(Exp.CrystalSymmetry)
    Exp.CrystalSymmetry = []; % no site splitting
  end
end
Exp.PowderSimulation = PowderSimulation;

% Partial ordering
if ~isempty(Exp.Ordering)
  if isnumeric(Exp.Ordering) && (numel(Exp.Ordering)==1) && isreal(Exp.Ordering)
    UserSuppliedOrderingFcn = 0;
    logmsg(1,'  partial order (built-in function, lambda = %g)',Exp.Ordering);
  elseif isa(Exp.Ordering,'function_handle')
    UserSuppliedOrderingFcn = 1;
    logmsg(1,'  partial order (user-supplied function)');
  else
    error('Exp.Ordering must be a single number or a function handle.');
  end
  if any(Sys.gStrain) || any(Sys.AStrain) || any(Sys.DStrain) || any(Sys.HStrain)
    error('Exp.Ordering and g/A/D/H strains cannot be used simultaneously.');
  end
end

% Temperature and non-equilibrium populations
NonEquiPops = 0;
if isfinite(Exp.Temperature)
  if numel(Exp.Temperature)==1
    msg = sprintf('  temperature %g K',Exp.Temperature);
  else
    msg = '  user-specified non-equilibrium populations';
    NonEquiPops = 1;
    if max(Exp.Temperature)==min(Exp.Temperature)
      error('Populations in Exp.Temperature cannot be all equal!');
    end
  end
else
  msg = '  no temperature';
end
logmsg(1,msg);
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
ObsoleteOptions = {'Convolution','Width'};
for k = 1:numel(ObsoleteOptions)
  if isfield(Opt,ObsoleteOptions{k}),
    error('Options.%s is obsolete. Please remove from code!',ObsoleteOptions{k});
  end
end

if isfield(Opt,'nSpline')
  error('Options.nSpline is obsolete. Use a second number in Options.nKnots instead, e.g. Options.nKnots = [19 5] for Options.nSpline = 5.');
end
if isfield(Opt,'Perturb')
  error('Options.Perturb is obsolete. Use Opt.Method=''perturb'' or Opt.Method=''hybrid'' instead.');
end

% Documented fields, pepper
DefaultOpt.Verbosity = 0;
DefaultOpt.Symmetry = 'auto';
DefaultOpt.SymmFrame = [];
DefaultOpt.Output = 'summed';
DefaultOpt.Method = 'matrix'; % 'matrix', 'eig', 'perturb1', 'perturb2'='perturb' 

% Undocumented fields, pepper
%DefaultOpt.nTRKnots = 3; % resfields
%DefaultOpt.Weak = []; % resfields
DefaultOpt.Smoothing = 2;
DefaultOpt.nKnotsMinimum = 10;
DefaultOpt.Debug = 0;
DefaultOpt.Intensity = 'on';
DefaultOpt.BruteForce = 0;
DefaultOpt.DirectAccumulation = 0;
DefaultOpt.PaddingMultiplier = 3; % for padding before convolution
DefaultOpt.ThetaRange = [];

nKnotsMatrix = [19 4];
nKnotsPerturb = [19 4];

Opt = adddefaults(Opt,DefaultOpt);

[Method,err] = parseoption(Opt,'Method',{'eig','matrix','perturb','perturb1','perturb2','hybrid'});
error(err);

UseEigenFields = (Method==1);
UsePerturbationTheory = (Method==3) || (Method==4) || (Method==5);
if (Opt.DirectAccumulation && ~UsePerturbationTheory)
  error('Opt.DirectAccumulation works only with perturbation theory.');
end

if ~isfield(Opt,'nKnots')
  if UsePerturbationTheory
    Opt.nKnots = nKnotsPerturb;
  else
    Opt.nKnots = nKnotsMatrix;
  end
end

if Opt.nKnots(1)<Opt.nKnotsMinimum
  error('Options.nKnots must not be less than %d. Please adjust!',Opt.nKnotsMinimum);
end
if numel(Opt.nKnots)<2
  if UsePerturbationTheory
    Opt.nKnots(2) = nKnotsPerturb(2);
  else
    Opt.nKnots(2) = nKnotsMatrix(2);
  end
end

% Parse string options.
[Opt.Output,err] = parseoption(Opt,'Output',{'summed','separate'});
error(err);
SummedOutput = (Opt.Output==1);

if isfield(Opt,'LineShape')
  disp('Options.LineShape is obsolete. Use System.lw or System.lwpp instead.');
  disp('- to specify a Gaussian broadening: System.lw = x.');
  disp('- to specify a Lorentzian broadening:   System.lw = [0 x].');
  disp('- to specify a Voigtian (Gaussian plus Lorentzian) broadening: System.lw = [x y].');
  error('Options.LineShape is obsolete. Use System.lw or System.lwpp instead.');
end

AnisotropicIntensities = parseoption(Opt,'Intensity',{'off','on'}) - 1;
Opt.Intensity = AnisotropicIntensities;

if strcmp(Opt.Symmetry,'auto'),
  Opt.Symmetry = [];
end

% Unclear which of the following versions is implemented in cw EPR spectrometers.
%Exp.deltaX = diff(Exp.Range)/Exp.nPoints; % excludes last value
Exp.deltaX = diff(Exp.Range)/(Exp.nPoints-1); % includes last value
xAxis = Exp.Range(1) + Exp.deltaX*(0:Exp.nPoints-1);


p_symandgrid;


%=======================================================================
%=======================================================================
%                   PEAK DATA COMPUTATION
%=======================================================================
%=======================================================================

logmsg(1,'-resonances--------------------------------------------');
MethodMsg{1} = 'eigenfields (Liouville space)';
MethodMsg{2} = 'adaptive segmentation (state space)';
MethodMsg{3} = 'second-order perturbation theory';
MethodMsg{4} = 'first-order perturbation theory';
MethodMsg{5} = 'second-order perturbation theory';
MethodMsg{6} = 'hybrid (matrix diagonalization for electron spin, perturbation for nuclei)';
logmsg(1,'  method: %s',MethodMsg{Method});

if (Method==1)
  
  % Eigenfield equation
  %------------------------------------------------------------------------
  AnisotropicWidths = 0;
  if StrainWidths
    logmsg(-inf,'WARNING: Options.Method: eigenfields method -> strains are ignored!');
    StrainWidths = 0;
  end

  Transitions = NaN;
  Exp1 = Exp;
  Exp1.Range = [0 1e8];
  Exp1.Orientations = [phi;theta;chi];
  
  logmsg(2,'  -entering eigfields----------------------------------');
  [Pdat,Idat] = eigfields(Sys,Exp1,Opt);
  logmsg(2,'  -exiting eigfields-----------------------------------');
  Wdat = [];
  Gdat = [];
  Transitions = [];
  
  if (nOrientations==1)
    Pdat = {Pdat};
    Idat = {Idat};
  end
  nReson = 0;
  for k = 1:nOrientations,
    nReson = nReson + numel(Pdat{k});
  end
  logmsg(1,'  %d resonance in total (%g per orientation)',nReson,nReson/nOrientations);
    
else
  
  % Matrix diagonalization and perturbation methods
  %------------------------------------------------------------------------
  
  if (~PowderSimulation)
    %if ~isfield(Opt,'Perturb'), Opt.Perturb = 0; end
  end
  
  % Set search range larger than requested field range
  Exp1 = Exp;
  Exp1.SearchRange = Exp1.Range + 0.2*diff(Exp.Range)*[-1 1];
  Exp1.SearchRange(Exp1.SearchRange<0) = 0;
 
  Exp1.Orientations = [phi;theta;chi];
  Exp1.AccumWeights = Weights;

  Opt.peppercall = 1;
  
  logmsg(2,'  -entering resfields*----------------------------------');
  switch Method
    case {2,6}
      [Pdat,Idat,Wdat,Transitions,Gdat] = resfields(Sys,Exp1,Opt);
    case {3,5}
      Opt.PerturbOrder = 2;
      [Pdat,Idat,Wdat,Transitions,spec] = resfields_perturb(Sys,Exp1,Opt);
    case 4
      Opt.PerturbOrder = 1;
      [Pdat,Idat,Wdat,Transitions,spec] = resfields_perturb(Sys,Exp1,Opt);
  end
  logmsg(2,'  -exiting resfields*-----------------------------------');
  
  nTransitions = size(Transitions,1);
  
  if isempty(Wdat)
    AnisotropicWidths = 0;
  else
    AnisotropicWidths = (max(Wdat(:))>0);
  end
  
end

if (~AnisotropicIntensities)
  logmsg(1,'  neglecting amplitude anisotropy, taking amplitude average');
  if ~UseEigenFields, Idat = mean(Idat(:))*ones(size(Idat)); end
end

if (~PowderSimulation) && ~isempty(Exp.CrystalSymmetry)
  nNewOrientations = size(Pdat,2);
  nSites = nNewOrientations/nOrientations;
else
  nSites = 1;
end

%=======================================================================
%=======================================================================


% Analysis of the resonance data
%------------------------------------------------------------------------
% If looping fields or out-of-range fields are encountered,
% skip interpolation and projection.
if ~UseEigenFields
  NaN_in_Pdat = any(isnan(Pdat),2);
else
  NaN_in_Pdat = 0;
end

if (Method~=6)
  LoopingTransitionsPresent = size(unique(Transitions,'rows'),1)<size(Transitions,1);
  if LoopingTransitionsPresent && PowderSimulation
    logmsg(0,'** Looping transitions found. Artifacts at coalescence points possible.');
  end
else
  % hybrid method: Transitions contains replicas of core sys transitions
  LoopingTransitionsPresent = 0;
end

%=======================================================================
%=======================================================================
%               INTERPOLATION AND SPECTRUM CONSTRUCTION
%=======================================================================
%=======================================================================
% The position/amplitude/width data from above are
%     (1) interpolated to get more of them, and then
%     (2) projected to obtain the spectrum
% The interpolation and projection algorithms depend
% on the symmetry of the data.

logmsg(1,'-absorption spectrum construction----------------------');

BruteForceSum = UseEigenFields | Opt.BruteForce;

if Opt.DirectAccumulation

  % bypass interpolation
  
elseif (~BruteForceSum)
  
  % Determine methods: projection/summation, interpolation on/off
  %-----------------------------------------------------------------------
  DoProjection = (~AnisotropicWidths) & (nOctants>=0);
  
  DoInterpolation = (Opt.nKnots(2)>1) & (nOctants>=0);
  
  % Preparations for projection
  %-----------------------------------------------------------------------
  if (DoProjection)
    msg = 'triangle/segment projection';
  else
    msg = 'summation';
    xT = 5e4;
    wT = xT/2.5; %<1e-8 at borders for Harmonic = -1
    Template = gaussian(0:2*xT-1,xT,wT,-1);
  end
  Text = {'single-crystal','isotropic','axial','nonaxial D2h','nonaxial C2h','','nonaxial Ci'};
  logmsg(1,'  %s, %s case',msg,Text{nOctants+3});
  
  % Preparations for interpolation
  %-----------------------------------------------------------------------
  if (DoInterpolation)
    % Set an option for the sparse tridiagonal matrix \ solver in global cubic
    % spline interpolation. This function needs some time, so it was taken
    % out of Matlab's original spline() function, which is called many times.
    spparms('autommd',0);
    % Interpolation parameters. 1st char: g global, l linear. 2nd char: order.
    if (nOctants==0), % axial symmetry: 1D interpolation
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
    msg = sprintf('  interpolation (factor %d, method %s/%s/%s)',Opt.nKnots(2),InterpMode{1},InterpMode{2},InterpMode{3});
  else
    Opt.nKnots(2) = 1;
    msg = '  interpolation off';
  end
  nfKnots = (Opt.nKnots(1)-1)*Opt.nKnots(2) + 1;
  logmsg(1,msg);
  
  % Pre-allocation of spectral array.
  %-----------------------------------------------------------------------
  if (SummedOutput)
    nRows = 1;
    msg = 'summed';
  else
    if (~PowderSimulation)
      nRows = nOrientations;
    else
      nRows = nTransitions;
    end
    msg = 'separate';
  end
  spec = zeros(nRows,Exp.nPoints);
  logmsg(1,'  spectrum array size: %dx%d (%s)',size(spec,1),size(spec,2),msg);
  
  
  %  DefaultOpt.SmoothingAlpha = 1.2;
  %  if (nOctants>=0)
  %    msg = [msg ' with gradient smoothing'];
  %    % Gradient-weighted line broadening for smoothing out simulation noise
  %    % (XSophe: mosaic misorientation model, Weihe: the best invention since sliced bread)
  %    nKnots(2) = max(Opt.nKnots(2),1);
  %    SmoothingPrefactor = (pi/2)/((Opt.nKnots(1)-1)*nKnots(2))*Opt.SmoothingAlpha;
  %    if isempty(Wdat), Wdat = 0; end
  %    Wdat = sqrt(Wdat.^2 + Gdat.^2*SmoothingPrefactor^2);
  %    AnisotropicWidths = 1;
  %  end
  
  
  % Spectrum construction
  %-----------------------------------------------------------------------
  if (~PowderSimulation)
    %=======================================================================
    % Single-crystal spectra
    %=======================================================================
    
    if (~AnisotropicIntensities), thisInt = ones(nTransitions,1); end
    if (~AnisotropicWidths), thisWid = zeros(nTransitions,1); end
    
    idx = 1;
    for iOri = 1:nOrientations
      for iSite = 1:nSites
        %logmsg(3,'  orientation %d of %d, site %d of %d',iOri,nOrientations,iSite,nSites);
        thisPos = Pdat(:,idx);
        if (AnisotropicIntensities), thisInt = Idat(:,idx); end
        if (AnisotropicWidths), thisWid = Wdat(:,idx); end
        
        thisspec = lisum1i(Template,xT,wT,thisPos,thisInt,thisWid,xAxis);
        
        if (SummedOutput)
          spec = spec + thisspec;
        else
          spec(iOri,:) = spec(iOri,:) + thisspec;
        end
        
        idx = idx + 1;
      end
    end
    
  elseif (nOctants==-1)
    
    %=======================================================================
    % Isotropic powder spectra
    %=======================================================================
    
    if (~AnisotropicIntensities), thisInt = 1; end
    if (~AnisotropicWidths), thisWid = 0; end
    
    for iTrans = 1:nTransitions
      %logmsg(3,'  transition %d of %d',iTrans,nTransitions);
      
      thisPos = Pdat(iTrans,:);
      if (AnisotropicIntensities), thisInt = Idat(iTrans,:); end
      if (AnisotropicWidths), thisWid = Wdat(iTrans,:); end
      
      thisspec = 4*pi*lisum1i(Template,xT,wT,thisPos,thisInt,thisWid,xAxis);
      
      if (SummedOutput)
        spec = spec + thisspec;
      else
        spec(iTrans,:) = thisspec;
      end
      
    end
    
  else
    
    %=======================================================================
    % Powder spectra: interpolation and accumulation/projection
    %=======================================================================
    Axial = (nOctants==0);
    if (Axial)
      if (DoInterpolation)
        [fphi,fthe] = sphgrid(Opt.Symmetry,nfKnots,'f');
      else
        fthe = theta;
      end
      fSegWeights = -diff(cos(fthe))*4*pi; % sum is 4*pi
      if ~isempty(Exp.Ordering)
        centreTheta = (fthe(1:end-1)+fthe(2:end))/2;
        if (UserSuppliedOrderingFcn)
          OrderingWeights = feval(Exp.Ordering,zeros(1,numel(centreTheta)),centreTheta);
          %OrderingWeights = Exp.Ordering(zeros(1,numel(centreTheta)),centreTheta);
          if any(OrderingWeights)<0, error('User-supplied orientation distribution gives negative values!'); end
          if max(OrderingWeights)==0, error('User-supplied orientation distribution is all-zero.'); end
        else
          OrderingWeights = exp(Exp.Ordering*plegendre(2,0,cos(centreTheta)));
        end
        fSegWeights = fSegWeights(:).*OrderingWeights(:);
        fSegWeights = 4*pi/sum(fSegWeights)*fSegWeights;
      elseif ~isempty(Opt.ThetaRange)
        centreTheta = (fthe(1:end-1)+fthe(2:end))/2;
        idx = (centreTheta<Opt.ThetaRange(1)) | (centreTheta>Opt.ThetaRange(2));
        fSegWeights(idx) = 0;
      end
      logmsg(1,'  total %d segments, %d transitions',numel(fthe)-1,nTransitions);
      
    else % nonaxial symmetry
      if (DoInterpolation)
        [fphi,fthe] = sphgrid(Opt.Symmetry,nfKnots,'f');
      else
        fthe = theta;
        fphi = phi;
      end
      [idxTri,Areas] = triangles(nOctants,nfKnots,ang2vec(fphi,fthe));
      if ~isempty(Exp.Ordering)
        centreTheta = mean(fthe(idxTri));
        if (UserSuppliedOrderingFcn)
          centrePhi = mean(fphi(idxTri));
          OrderingWeights = feval(Exp.Ordering,centrePhi,centreTheta);
          %OrderingWeights = Exp.Ordering(centrePhi,centreTheta);
          if any(OrderingWeights)<0, error('User-supplied orientation distribution gives negative values!'); end
          if max(OrderingWeights)==0, error('User-supplied orientation distribution is all-zero.'); end
        else
          OrderingWeights = exp(Exp.Ordering*plegendre(2,0,cos(centreTheta)));
        end
        Areas = Areas(:).*OrderingWeights(:);
        Areas = 4*pi/sum(Areas)*Areas;
      elseif ~isempty(Opt.ThetaRange)
        centreTheta = mean(fthe(idxTri));
        idx = (centreTheta<Opt.ThetaRange(1)) | (centreTheta>Opt.ThetaRange(2));
        Areas(idx) = 0;
      end
      logmsg(1,'  total %d triangles (%d orientations), %d transitions',size(idxTri,2),numel(fthe),nTransitions);
    end
    
    if (~AnisotropicIntensities), fInt = ones(size(fthe)); end
    if (~AnisotropicWidths), fWid = zeros(size(fthe)); end
    
    minBroadening = inf;
    nBroadenings = 0;
    sumBroadenings = 0;
    
    for iTrans = 1:nTransitions
      
      % Interpolation
      %------------------------------------------------------
      %LoopTransition = any(isnan(Pdat(iTrans,:)));
      LoopTransition = 0;
      InterpolateThis = (DoInterpolation & ~LoopTransition);
      if (InterpolateThis)
        fPos = esintpol(Pdat(iTrans,:),InterpParams,Opt.nKnots(2),InterpMode{1},fphi,fthe);
        if (AnisotropicIntensities)
          fInt = esintpol(Idat(iTrans,:),InterpParams,Opt.nKnots(2),InterpMode{2},fphi,fthe);
        end
        if (AnisotropicWidths)
          fWid = esintpol(Wdat(iTrans,:),InterpParams,Opt.nKnots(2),InterpMode{3},fphi,fthe);
        end
      else
        fPos = Pdat(iTrans,:);
        if (AnisotropicIntensities), fInt = Idat(iTrans,:); end
        if (AnisotropicWidths), fWid = Wdat(iTrans,:); end
      end
      
      msg1 = '';
      if (~NonEquiPops) && any(fInt(:)<0), msg1 = 'intensities'; end
      if any(fWid(:)<0), msg1 = 'widths'; end
      if ~isempty(msg1)
        error('Negative %s encountered! Please report!',msg1);
      end
      
      % Summation or projection
      %------------------------------------------------------
      if (DoProjection && ~LoopTransition)
        if (Axial)
          thisspec = projectzones(fPos,fInt,fSegWeights,xAxis);
        else
          thisspec = projecttriangles(idxTri,Areas,fPos,fInt,xAxis);
        end
        % minBroadening = ?
      else % do summation
        if (Axial)
          fPosC = (fPos(1:end-1) + fPos(2:end))/2;
          fIntC = fSegWeights.*(fInt(1:end-1) + fInt(2:end))/2;
          fSpread = abs(fPos(1:end-1) - fPos(2:end));
          fWidM  = (fWid(1:end-1) + fWid(2:end))/2;
          c1 = 1.57246; c2 = 18.6348;
        else
          fPosSorted = sort(fPos(idxTri),1);
          fPosC = mean(fPosSorted,1);
          fIntC = Areas.*mean(fInt(idxTri),1);
          fSpread = fPosSorted(3,:) - fPosSorted(1,:);
          fWidM = mean(fWid(idxTri),1);
          c1 = 2.8269; c2 = 42.6843;
        end
        Lambda = fWidM./fSpread;
        gam = 1./sqrt(c1*Lambda.^2 + c2*Lambda.^4);
        fWidC = fWidM.*(1 + Opt.Smoothing*gam);
        
        thisspec = lisum1i(Template,xT,wT,fPosC,fIntC,fWidC,xAxis);
        
        minBroadening = min(minBroadening,min(Lambda));
        sumBroadenings = sumBroadenings + sum(Lambda);
        nBroadenings = nBroadenings + numel(Lambda);
      end
      
      % Accumulate subspectra
      %----------------------------------------------------------
      if (SummedOutput)
        spec = spec + thisspec;
      else
        spec(iTrans,:) = thisspec;
      end
      
    end % for iTrans
    
    if (~DoProjection)
      logmsg(1,'  Smoothness: overall %0.4g, worst %0.4g\n   (<0.5: probably bad, 0.5-3: ok, >3: overdone)',sumBroadenings/nBroadenings,minBroadening);
    end
    
  end
  %=======================================================================
  
else % if ~BruteForceSum else ...
  
  logmsg(1,'  no interpolation',nOrientations);
  logmsg(1,'  constructing stick spectrum');
  logmsg(1,'  summation over %d orientations',nOrientations);
  spec = zeros(1,Exp.nPoints);
  prefactor = (Exp.nPoints-1)/(Exp.Range(2)-Exp.Range(1));
  for k = 1:nOrientations
    if iscell(Pdat)
      p = round(1+prefactor*(Pdat{k}-Exp.Range(1)));
      i = Idat{k};
    else
      p = round(1+prefactor*(Pdat(:,k)-Exp.Range(1)));
      i = Idat(:,k);
    end
    OutOfRange = (p<1) | (p>Exp.nPoints);
    p(OutOfRange) = [];
    i(OutOfRange) = [];
    if (AnisotropicIntensities)
      spec = spec + full(sparse(1,p,Weights(k)*i,1,Exp.nPoints));
    else
      spec = spec + full(sparse(1,p,Weights(k),1,Exp.nPoints));
    end
  end
  spec = spec/Exp.deltaX;
  
end
%=======================================================================





%=======================================================================
%                         Final activities
%=======================================================================
logmsg(1,'-final-------------------------------------------------');

% Combine branches of looping transitions if separate output
%-----------------------------------------------------------------------
if (PowderSimulation)
  if (~SummedOutput) && LoopingTransitionsPresent
    [Transitions,unused,idx] = unique(Transitions,'rows');
    nTransitions = size(Transitions,1);
    newspec = zeros(nTransitions,Exp.nPoints);
    for k = 1:length(idx)
      newspec(idx(k),:) = newspec(idx(k),:) + spec(k,:);
    end
    spec = newspec;
    clear newspec;
  end
end

% Convolution with line shape.
%-----------------------------------------------------------------------
if (ConvWidth)
  logmsg(1,'  harmonic %d: using convolution',Exp.Harmonic);
  lwG = Sys.lw(1);
  lwL = Sys.lw(2);
  Harmonic2Do = Exp.Harmonic;

  % Add padding to left and right of spectral range
  % to reduce convolution artifacts
  nPad = 0;
  exceedsLowerLimit = any(spec(:,1)~=0);
  exceedsHigherLimit = any(spec(:,end)~=0);
  if exceedsLowerLimit
    if exceedsHigherLimit
      logmsg(0,'** Spectrum exceeds field range. Artifacts at lower and upper field limits possible.');
    else
      logmsg(0,'** Spectrum exceeds field range. Artifacts at lower field limit possible.');
    end
  else
    if exceedsHigherLimit
      logmsg(0,'** Spectrum exceeds field range. Artifacts at upper field limit possible.');
    end
  end
  if  exceedsLowerLimit || exceedsHigherLimit
    nPad = round(max([lwG lwL])/Exp.deltaX*Opt.PaddingMultiplier);
    spec = [repmat(spec(:,1),1,nPad) spec]; % left padding
    spec = [spec repmat(spec(:,end),1,nPad)]; % right padding
  end

  % Gaussian broadening
  if (lwG>0)
    logmsg(1,'  convoluting with Gaussian, FWHM %g mT, derivative %d',lwG,Harmonic2Do);
    if min(size(spec))==1, fwhm = [lwG 0]; else fwhm = [0 lwG]; end
    spec = convspec(spec,Exp.deltaX,fwhm,Harmonic2Do,1);
    Harmonic2Do = 0;
  end

  % Lorentzian broadening
  if (lwL>0)
    logmsg(1,'  convoluting with Lorentzian, FWHM %g mT, derivative %d',lwL,Harmonic2Do);
    if min(size(spec))==1, fwhm = [lwL 0]; else fwhm = [0 lwL]; end
    spec = convspec(spec,Exp.deltaX,fwhm,Harmonic2Do,0,Exp.mwPhase);
    Harmonic2Do = 0;
  else
    if (Exp.mwPhase~=0)
      error('For Exp.mwPhase different from zero, please specify a Lorentzian broadening in Sys.lwpp or Sys.lw.');
    end
  end

  % Remove padding
  if (nPad>0)
    spec(:,1:nPad) = [];
    spec(:,Exp.nPoints+1:end) = [];
  end

else
  
  if (Exp.Harmonic>0)
    logmsg(1,'  harmonic %d: using differentiation',Exp.Harmonic);
    for h = 1:Exp.Harmonic
      dspec = diff(spec,[],2)/Exp.deltaX;
      spec = (dspec(:,[1 1:end]) + dspec(:,[1:end end]))/2;
    end
  else
    logmsg(1,'  harmonic 0: absorption spectrum');
  end

end

% Field modulation
%-----------------------------------------------------------------------
if (Exp.ModAmp>0)
  logmsg(1,'  applying field modulation');
  spec = fieldmod(xAxis,spec,Exp.ModAmp);
else
  % derivatives already included in convolutions etc.
end

% Assign output.
%-----------------------------------------------------------------------
switch (nargout)
  case 0,
  case 1, varargout = {spec};
  case 2, varargout = {xAxis,spec};
  case 3, varargout = {xAxis,spec,Transitions};
end

% Report performance.
%-----------------------------------------------------------------------
[Hours,Minutes,Seconds] = elapsedtime(StartTime,clock);
msg = sprintf('cpu time %dh%dm%0.3fs',Hours,Minutes,Seconds);
logmsg(0.5,msg);

logmsg(1,'=end=pepper=======%s=================\n',datestr(now));

clear global EasySpinLogLevel

return
