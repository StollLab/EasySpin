% salt  ENDOR spectra simulation
%
%   salt(Sys,Exp)
%   salt(Sys,Exp,Opt)
%   spec = salt(...)
%   [rf,spec] = salt(...)
%   [rf,spec,info] = salt(...)
%
%   Input:
%   - Sys: spin system structure specification
%       S, g, gFrame, Nucs, A, AFrame, Q, QFrame etc.
%       lwEndor     FWHM Endor line width [Gaussian,Lorentzian]
%   - Exp: experiment specification
%       mwFreq        spectrometer frequency, in GHz
%       Field         magnetic field, in mT
%       Range         radiofrequency range [low,high], in MHz
%       nPoints       number of points
%       Temperature   temperature of the sample, by default off (NaN)
%       ExciteWidth   ENDOR excitation width, FWHM, in MHz
%   - Opt: simulation options
%       Verbosity           0, 1, 2
%       Method        'matrix', 'perturb1', 'perturb2'='perturb'
%       separate      '', 'components', 'transitions', 'sites', 'orientations'
%       GridSize      grid size;  N1, [N1 Ninterp]
%       Transitions, Threshold, GridSymmetry
%       Intensity, Enhancement, Sites
%
%   Output:
%   - rf:     the radiofrequency axis, in MHz
%   - spec:   the spectrum or spectra
%   - info:    structure with details of the calculation
%
%   If no output argument is given, the simulated spectrum is plotted.

function varargout = salt(Sys,Exp,Opt)

if nargin==0, help(mfilename); return; end

% Check expiry date
error(eschecker);

% Check Matlab version.
warning(chkmlver);

% Get time for performance prompt at the end.
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
% is used in logmsg(), which does the display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel;
EasySpinLogLevel = Opt.Verbosity;

%==================================================================
% Loop over species and isotopologues
%==================================================================
SweepAutoRange = (~isfield(Exp,'Range') || isempty(Exp.Range)) && ...
  (~isfield(Exp,'CenterSweep') || isempty(Exp.CenterSweep));

if ~isfield(Opt,'IsoCutoff'), Opt.IsoCutoff = 1e-4; end

% Process Opt.separate
if ~isfield(Opt,'separate'), Opt.separate = ''; end
[separateOutput,err] = parseoption(Opt,'separate',{'','components','transitions','orientations','sites'});
error(err);
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
  [xAxis,spec,info] = compisoloop(@salt,Sys,Exp,Opt,SweepAutoRange,thirdOutput,separateComponentSpectra);
  
  % Output and plotting
  switch nargout
    case 0
      salt_plot(xAxis,spec,Exp);
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

logmsg(1,'=begin=salt=======%s=================',char(datetime));
logmsg(2,'  log level %d',EasySpinLogLevel);
logmsg(1,'-general-----------------------------------------------');


% Processing spin system structure
%==========================================================================

[Sys,err] = validatespinsys(Sys);
error(err);
if Sys.MO_present, error('salt does not support Sys.Ham* parameters.'); end
if any(Sys.L(:)), error('salt does not support Sys.L.'); end

if any(Sys.n>1)
  error('salt does not support sets of equivalent nuclei (Sys.n>1).');
end

ConvolutionBroadening = any(Sys.lwEndor>0);

logmsg(1,'  system with %d electron spin(s), %d nuclear spin(s), total %d states',...
  Sys.nElectrons,Sys.nNuclei,Sys.nStates);

if isfield(Sys,'ExciteWidth')
  error('You gave Sys.ExciteWidth, but it should be Exp.ExciteWidth. Please correct.');
end
%==========================================================================


% Method: 'matrix' or 'perturb', 'perturb1', 'perturb2'
if ~isfield(Opt,'Method'), Opt.Method = 'matrix'; end
Method = parseoption(Opt,'Method',{'matrix','perturb','perturb1','perturb2'});
switch Method
  case 1, Opt.PerturbOrder = 0; Method = 1;
  case 2, Opt.PerturbOrder = 1; Method = 2;
  case 3, Opt.PerturbOrder = 1; Method = 2;
  case 4, Opt.PerturbOrder = 2; Method = 2;
end

%=======================================================================
% Experiment structure, contains experimental settings
%=======================================================================

% Documented fields and their defaults (mandatory parameters are set to NaN)
DefaultExp.Field = NaN;
DefaultExp.mwFreq = NaN;
DefaultExp.Range = NaN;
DefaultExp.CenterSweep = NaN;
DefaultExp.nPoints = 1024;
DefaultExp.Temperature = NaN;
DefaultExp.ExciteWidth = Inf;

DefaultExp.Ordering = [];
DefaultExp.SampleRotation = [];
DefaultExp.SampleFrame = [];
DefaultExp.CrystalSymmetry = '';
DefaultExp.MolFrame = [];

% Undocumented and unused fields
DefaultExp.Harmonic = 0;

% Make user-given fields case-insensitive.
Exp = adddefaults(Exp,DefaultExp);

% Check for obsolete fields
if isfield(Exp,'Orientations')
  error('Exp.Orientations is no longer supported, use Exp.SampleFrame instead.');
end
if isfield(Exp,'CrystalOrientation')
  error('Exp.CrystalOrientation is no longer supported, use Exp.SampleFrame instead.');
end

% Photoselection is not supported
if isfield(Exp,'lightBeam') && ~isempty(Exp.lightBeam)
  error('Photoselection (via Exp.lightBeam) is not supported.')
end

% Error if vital parameters are missing.
if isnan(Exp.Field), error('Experiment.Field is missing!'); end
  
if isnan(Exp.mwFreq)
  if ~isinf(Exp.ExciteWidth)
    error('Experiment.mwFreq is missing! It is needed, since Exp.ExciteWidth is given.');
  end
end

AutoRange = isnan(Exp.Range) & isnan(Exp.CenterSweep);
if Method==1
%  if AutoRange, error('Cannot automatically determine rf range. Please specify Exp.Range or Exp.CenterSweep.'); end
end

% Check both CenterSweep and Range, prefer CenterSweep
if ~isnan(Exp.CenterSweep)
  if ~isnan(Exp.Range)
    logmsg(0,'Using Experiment.CenterSweep and ignoring Experiment.Range.');
  end
  if Exp.CenterSweep(2)<=0
    error('Sweep range in Exp.CenterSweep must be positive.');
  end
  Exp.Range = Exp.CenterSweep(1) + [-1 1]*Exp.CenterSweep(2)/2;
  if Exp.Range(1)<0
    error('Start value resulting from Exp.CenterSweep is negative (%g).',Exp.Range(1));
  end
end

% Check Exp.Range
if ~isnan(Exp.Range)
  if any(diff(Exp.Range)<=0) || ...
      any(~isfinite(Exp.Range)) || ...
      ~isreal(Exp.Range) || ...
      any(Exp.Range<0) || ...
      (numel(Exp.Range)~=2)
    error('Experiment.Range is not valid!');
  end
end

% Automatic on-resonance setting: (1) isotropic g, (2) mwFreq==0
if Exp.mwFreq==0 && (all(diff(Sys.g)==0))
  Exp.mwFreq = mean(Sys.g)*bmagn*Exp.Field*1e-3/planck/1e9;
end

% Detect sample type (disordered, partially ordered, crystal)
[Exp,Opt] = p_sampletype(Exp,Opt);
crystalSample = Opt.crystalSample;

if Opt.partiallyOrderedSample
  error('Partially ordered samples are not implemented in salt.')
end

if EasySpinLogLevel>=1
  if ~AutoRange
    msg = sprintf('field %g mT, rf range [%g %g] MHz, %d points',...
      Exp.Field,Exp.Range(1),Exp.Range(2),Exp.nPoints);
  else
    msg = sprintf('field %g mT, automatic rf range, %d points',...
      Exp.Field,Exp.nPoints);
  end
  if ~isnan(Exp.mwFreq)
    msg = sprintf('mw %g GHz, %s',Exp.mwFreq,msg);
  else
    msg = sprintf('mw not given, %s',msg);
  end
  logmsg(1,['  ' msg]);
end

nonEquiPops = false;
if isfinite(Exp.Temperature)
  msg = sprintf('  temperature %g K',Exp.Temperature);
else
  msg = '  no temperature';
end
logmsg(1,msg);

if Exp.Harmonic>0 && all(Sys.lwEndor==0)
  logmsg(-Inf,'WARNING: Cannot compute derivative of a stick spectrum! Returning absorption spectrum.');
  logmsg(-Inf,'WARNING:   Add a line width (strain or convolution) or set Harmonic=0.');
  Exp.Harmonic = 0;
end

if isfield(Exp,'HStrain')
  error('You gave Exp.HStrain, but it should be Sys.HStrain. Please correct.');
end

%==========================================================================


% Processing Options structure
%==========================================================================
% Documented fields, endorfrq (for defaults see there!)
%DefaultOpt.Transitions = []; % endorfrq
%DefaultOpt.Threshold = 5e-4; % endorfrq
%DefaultOpt.Nuclei = []; % endorfrq
%DefaultOpt.Intensity = 'on'; % endorfrq
%DefaultOpt.Enhancement = 'off'; % endorfrq
%DefaultOpt.Sites = []; % endorfrq
% Documented fields, salt
if Method==2
  DefaultOpt.GridSize = [61 0];
else
  DefaultOpt.GridSize = [31 3];
end
DefaultOpt.GridSymmetry = '';
DefaultOpt.GridFrame = [];
DefaultOpt.separate = '';

% Undocumented fields
DefaultOpt.Smoothing = 2;
DefaultOpt.minEffKnots = 5;

% Undocumented fields
DefaultOpt.OriWeights = [];
DefaultOpt.OriThreshold = 1e-4;
DefaultOpt.IncludeOriWeights = 1;
DefaultOpt.OriPreSelect = [];
DefaultOpt.PaddingMultiplier = 3; % for padding before convolution

% Supplement the user-supplied structure with the defaults
Opt = adddefaults(Opt,DefaultOpt);


if Opt.GridSize(1)<2
  error('Opt.GridSize(1) must be 2 or larger.');
end
if numel(Opt.GridSize)<2
  Opt.GridSize(2) = DefaultOpt.GridSize(2);
end

% Some compatibility checks for separate spectra output (Opt.separate)
if ~crystalSample
  if separateOrientationSpectra
    error(sprintf('\nCannot return separate orientations for non-crystal spectra (Opt.separate=''orientations'').\nUse other setting for Opt.separate.\n'));
  end
  if separateSiteSpectra
    error(sprintf('\nCannot return separate sites for non-crystal spectra (Opt.separate=''sites'').\nUse other setting for Opt.separate.\n'));
  end
else
  if separateTransitionSpectra
    error(sprintf('\n  Cannot return separate transitions for crystal spectra (Opt.separate=''transitions'').\n  Use other setting for Opt.separate.\n'));
  end
end


suppliedOriWeights = ~isempty(Opt.OriWeights);
if suppliedOriWeights
  Opt.OriWeights = Opt.OriWeights(:).';
end
if any(Opt.OriPreSelect) && suppliedOriWeights
  error('You cannot use OriPreSelect and supply OriWeights simultaneously!');
end

% Error if obsolete option is specified
if isfield(Opt,'nKnots')
  error('Options.nKnots is obsolete. Use Options.GridSize instead, e.g. Options.GridSize = 91.');
end

if isfield(Opt,'Symmetry')
  error('Options.Symmetry is obsolete. Use Options.GridSymmetry instead, e.g. Options.GridSymmetry = ''D2h''.');
end

if isfield(Opt,'Convolution')
  error('Options.Convolution is obsolete! Please remove from code!');
end

if isfield(Opt,'Width')
  error('Options.Width is obsolete! Please remove from code.');
end

if isfield(Opt,'nSpline')
  error('Options.nSpline is obsolete. Use a second number in Options.GridSize instead, e.g. Options.GridSize = [19 5] for Options.nSpline = 5.');
end

if isfield(Opt,'LineShape')
  error('Options.LineShape is obsolete. Use System.lwEndor instead.');
end

if isfield(Opt,'Output')
  error('Options.Output is obsolete. Use Options.separate instead.');
end

%==========================================================================

%====================================================================
% Special case: no nuclei
%--------------------------------------------------------------
if Sys.nNuclei==0
  Exp.deltaX = diff(Exp.Range)/(Exp.nPoints-1);
  xAxis = Exp.Range(1) + (0:Exp.nPoints-1)*Exp.deltaX;
  spec = zeros(1,Exp.nPoints);
  Transitions = [];
  switch (nargout)
    case 0
      cla
      plot(xAxis,spec);
      axis tight
      xlabel('frequency (MHz)');
      ylabel('intensity (arb.u.)');
      if ~isnan(Exp.mwFreq)
        title(sprintf('ENDOR at %0.8g GHz and %0.8g mT',Exp.mwFreq,Exp.Field));
      else
        title(sprintf('ENDOR at %0.8g mT',Exp.Field));
      end
    case 1
      varargout = {spec};
    case 2
      varargout = {xAxis,spec};
    case 3
      info.Transitions = Transitions;
      varargout = {xAxis,spec,info};
  end
  return
end
%====================================================================

% Set up orientational grid
[Exp,Opt,nOrientations] = p_gridsetup(Sys,Exp,Opt);


%==========================================================================
% Orientation pre-selection
%==========================================================================
% Based on a reduced spin system, weights are computed for orientation.
% They indicate whether any EPR transition is excited for that orientation
% or not. Based on the weight for a particular orientation and on
% Opt.OriThreshold, endorfrq will decide whether to compute the ENDOR
% spectrum.

Opt.DoPreSelection = ~crystalSample & Opt.OriPreSelect & ...
  Opt.OriThreshold>0 & isempty(Opt.OriWeights);

if Opt.DoPreSelection
  logmsg(1,'  computing orientations which fall into excitation window');
  OriSys = anisosubsys(Sys);
  OriExp = Exp;
  OriExp.Orientations = [Exp.phi;Exp.theta];
  Opt.OriWeights = orisel(OriSys,OriExp);
  Opt.OriWeights = Opt.OriWeights/max(Opt.OriWeights);
  str = '  %g percent of all orientations above relative threshold %g';
  logmsg(1,sprintf(str,100*sum(Opt.OriWeights>Opt.OriThreshold)/nOrientations,Opt.OriThreshold));
  %plot(Exp.theta*180/pi,Opt.OriWeights)
end
%==========================================================================



%=======================================================================
%=======================================================================
%                   PEAK DATA COMPUTATION
%=======================================================================
%=======================================================================

logmsg(1,'-resonances--------------------------------------------');
MethodName = {'  method: matrix diagonalization','  method: perturbation theory'};
logmsg(1,MethodName{Method});
logmsg(2,'  -endorfrq start-----------------------------------');
switch Method
  case 1, [Pdat,Idat,Transitions,Info] = endorfrq(Sys,Exp,Opt);
  case 2, [Pdat,Idat,Transitions,Info] = endorfrq_perturb(Sys,Exp,Opt);
end
logmsg(2,'  -endorfrq end-------------------------------------');
nTransitions = size(Pdat,1);
anisotropicWidths = false;

if nTransitions==0
  error(['No ENDOR resonances between %g and %g MHz (%g mT).\n'...
      'Check frequency range, magnetic field and transition selection.'],...
    Exp.Range(1),Exp.Range(2),Exp.Field);
end
logmsg(1,'  %d transitions with resonances in range',nTransitions);
logmsg(2,'  positions min %g MHz, max %g MHz',min(Pdat(:)),max(Pdat(:)));

anisotropicIntensities = numel(Idat)>1;
if anisotropicIntensities
  logmsg(2,'  amplitudes min %g, max %g',min(Idat(:)),max(Idat(:)));
end

if suppliedOriWeights
  logmsg(1,'  user supplied OriWeights: including them in intensity');
  if Opt.IncludeOriWeights
    for iT = 1:nTransitions
      Idat(iT,:) = Idat(iT,:).*Opt.OriWeights;
    end
  end
end

if crystalSample && ~isempty(Exp.CrystalSymmetry)
  nSites = numel(Pdat)/nTransitions/nOrientations;
else
  nSites = 1;
end

if AutoRange
  Range = [min(Pdat(:)) max(Pdat(:))]+[-1 1]*max(Sys.lwEndor)*2;
  Range = Range + diff(Range)*[-1 1]*0.2;
  Range(1) = max(Range(1),0); 
  Exp.Range = Range;
end
Exp.deltaX = diff(Exp.Range)/(Exp.nPoints-1); % includes last value
xAxis = Exp.Range(1) + (0:Exp.nPoints-1)*Exp.deltaX;

%=======================================================================
%=======================================================================



% Guard against strong orientation selectivity.
%----------------------------------------------------------------------------
% If the width of the EPR spectrum is much greater than the excitation width,
% interpolation and projection is not advisable.

if Info.Selectivity>0
  logmsg(1,'  orientation selection: %g (<1 very weak, 1 weak, 10 strong, >10 very strong)',Info.Selectivity);
  GridTooCoarse = (Opt.GridSize(1)/Opt.minEffKnots<Info.Selectivity);
  if GridTooCoarse && ~crystalSample
    fprintf('  ** Warning: Strong orientation selection ********************************\n');
    fprintf('  Only %0.1f orientations within excitation window.\n',Opt.GridSize(1)/Info.Selectivity);
    fprintf('  Spectrum might be inaccurate!\n');
    fprintf('  To remedy, do one (or both) of the following:\n');
    fprintf('  - Increase Opt.GridSize (currently %d)\n',Opt.GridSize(1));
    fprintf('  - Increase Exp.ExciteWidth (currently %g MHz)\n',Exp.ExciteWidth);
    fprintf('  *************************************************************************\n');
  end
else
  logmsg(1,'  orientation selection: none');
end

%=======================================================================
%=======================================================================
%    SPECTRUM CONSTRUCTION (incl. INTERPOLATION and PROJECTION)
%=======================================================================
%=======================================================================
% The position/amplitude/width data from above are
%     (1) interpolated to get more of them, and then
%     (2) projected to obtain the spectrum
% The interpolation and projection algorithms depend
% on the symmetry of the data.

logmsg(1,'-absorption spectrum construction----------------------');

NaN_in_Pdat = any(isnan(Pdat),2);
if any(NaN_in_Pdat)
  Opt.GridSize(2) = 1;
end

BruteForceSum = false;
axialGrid = Opt.nOctants==0;
usingGrid = Opt.nOctants>=0;

if ~BruteForceSum

  % Preparations for summation/projection
  %-----------------------------------------------------------------------
  doProjection = usingGrid && ~anisotropicWidths;
  if doProjection
    msg = 'triangle/segment projection';
  else
    msg = 'summation';
    Template.x = 5e4;
    Template.lw = Template.x/2.5; %<1e-8 at borders for Harmonic = -1
    Template.y = gaussian(0:2*Template.x-1,Template.x,Template.lw,-1);
  end
  logmsg(1,'  %s',msg);

  % Preparations for interpolation
  %-----------------------------------------------------------------------
  doInterpolation = usingGrid && Opt.GridSize(2)>1;
  if doInterpolation
    % Set an option for the sparse tridiagonal matrix \ solver in global cubic
    % spline interpolation. This function needs some time, so it was taken
    % out of Matlab's original spline() function, which is called many times.
    spparms('autommd',0);
    % Interpolation parameters. 1st char: g global, l linear. 2nd char: order.
    if axialGrid  % axial symmetry: 1D interpolation
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
  if crystalSample
    %=======================================================================
    % Single-crystal spectra
    %=======================================================================

    if ~anisotropicIntensities, thisInt = ones(nTransitions,1); end
    if ~anisotropicWidths, thisWid = zeros(nTransitions,1); end

    iOriSite = 1;  % index into Pdat/Idat/Wdat
    spcidx = 0;  % index into spectral output array
    for iOri = 1:nOrientations
      if separateSiteSpectra, spcidx = 0;
      elseif separateOrientationSpectra, spcidx = spcidx+1;
      else, spcidx = 1;
      end
      for iSite = 1:nSites
        if separateSiteSpectra, spcidx = spcidx + 1; end
        %logmsg(3,'  orientation %d of %d',iOri,nOrientations);

        thisPos = Pdat(:,iOriSite);
        if anisotropicIntensities, thisInt = Idat(:,iOriSite); end
        %if anisotropicWidths, thisWid = Wdat(:,iOriSite); end

        thisspec = lisum1i(Template.y,Template.x,Template.lw,thisPos,thisInt,thisWid,xAxis);
        thisspec = thisspec/nSites/nOrientations;
        thisspec = (2*pi)*thisspec; % for consistency with powder spectra (integral of chi)
        thisspec = Exp.OriWeights(iOri)*thisspec; % for consistency with powder spectra (integral over phi,theta)

        if ~separateSiteSpectra && ~separateOrientationSpectra
          spec = spec + thisspec;
        else
          spec(spcidx,:) = spec(spcidx,:) + thisspec;
        end

        iOriSite = iOriSite + 1;
      end
    end

  elseif ~usingGrid

    %=======================================================================
    % Isotropic disordered spectra
    %=======================================================================

    if ~anisotropicIntensities, thisInt = Idat; end
    if ~anisotropicWidths, thisWid = 0; end

    spcidx = 0;
    for iTrans = 1:nTransitions
      %logmsg(3,'  transition %d of %d',iTrans,nTransitions);

      thisPos = Pdat(iTrans,:);
      if anisotropicIntensities, thisInt = Idat(iTrans,:); end
      %if AnisotropicWidths, thisWid = Wdat(iTrans,:); end

      thisspec = lisum1i(Template.y,Template.x,Template.lw,thisPos,thisInt,thisWid,xAxis);
      thisspec = thisspec*(2*pi); % orientation integral (chi)
      thisspec = Exp.OriWeights*thisspec; % orientation integral (phi,theta)

      if ~separateTransitionSpectra
        spec = spec + thisspec;
      else
        spcidx = spcidx + 1;
        spec(spcidx,:) = thisspec;
      end

    end

  else

    %=======================================================================
    % Disordered spectra: interpolation and accumulation/projection
    %=======================================================================
    if axialGrid
      if doInterpolation
        grid = sphgrid(Opt.GridSymmetry,nfKnots);
        fthe = grid.theta;
        fphi = grid.phi;
      else
        fthe = Exp.theta;
      end
      fSegWeights = -diff(cos(fthe))*4*pi; % sum is 4*pi

      logmsg(1,'  total %d segments, %d transitions',numel(fthe)-1,nTransitions);

    else % nonaxial symmetry
      if doInterpolation
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

      logmsg(1,'  total %d triangles (%d orientations), %d transitions',size(idxTri,2),numel(fthe),nTransitions);
    end

    if ~anisotropicIntensities, fInt = ones(size(fthe)); end
    if ~anisotropicWidths, fWid = zeros(size(fthe)); end

    minBroadening = Inf;
    nBroadenings = 0;
    sumBroadenings = 0;
    spcidx = 0;

    for iTrans = 1:nTransitions

      % Interpolation
      %------------------------------------------------------
      if doInterpolation
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
      if ~nonEquiPops && any(fInt(:)<0)
        % Try to correct underswings due to (very) small numerical errors
        maxInt = max(fInt);
        fInt((fInt<0)&(fInt/maxInt>-1e-6)) = 0;
        if any(fInt(:)<0)
          msg1 = 'intensities';
        end
      end
      if any(fWid(:)<0), msg1 = 'widths'; end
      if ~isempty(msg1)
        error('Negative %s encountered! Please report!',msg1);
      end

      % Summation or projection
      %------------------------------------------------------
      if doProjection
        if axialGrid
          thisspec = projectzones(fPos,fInt,fSegWeights,xAxis);
        else
          thisspec = projecttriangles(idxTri,Areas,fPos,fInt,xAxis);
        end
        % minBroadening = ?
      else
        if axialGrid
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

        thisspec = lisum1i(Template.y,Template.x,Template.lw,fPosC,fIntC,fWidC,xAxis);

        minBroadening = min(minBroadening,min(Lambda));
        sumBroadenings = sumBroadenings + sum(Lambda);
        nBroadenings = nBroadenings + numel(Lambda);
      end

      thisspec = thisspec*(2*pi); % orientation integal over chi (0..2*pi)

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

else % if ~BruteForceSum else ...

  logmsg(1,'  no interpolation',nOrientations);
  logmsg(1,'  constructing stick spectrum');
  logmsg(1,'  summation over %d orientations',nOrientations);
  spec = zeros(1,Exp.nPoints);
  prefactor = (Exp.nPoints-1)/(Exp.Range(2)-Exp.Range(1));
  for k = 1:nOrientations
    p = round(1+prefactor*(Pdat{k}-Exp.Range(1)));
    OutOfRange = (p<1) | (p>Exp.nPoints);
    p(OutOfRange) = [];
    Idat{k}(OutOfRange) = [];
    if anisotropicIntensities
      spec = spec + full(sparse(1,p,Exp.OriWeights(k)*Idat{k},1,Exp.nPoints));
    else
      spec = spec + full(sparse(1,p,Exp.OriWeights(k),1,Exp.nPoints));
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

% Convolution with line shape.
%-----------------------------------------------------------------------
logmsg(1,'  harmonic %d',Exp.Harmonic);
if ConvolutionBroadening
  fwhmG = Sys.lwEndor(1);
  fwhmL = Sys.lwEndor(2);
  if fwhmL>0
    HarmonicL = Exp.Harmonic;
    HarmonicG = 0;
  else
    HarmonicL = 0;
    HarmonicG = Exp.Harmonic;
  end

  % Add padding to left and right of spectral range
  % to reduce convolution artifacts
  RangePadding = true;
  if RangePadding
    exceedsLowerLimit = any(spec(:,1)~=0);
    exceedsHigherLimit = any(spec(:,end)~=0);
    if exceedsLowerLimit
      if exceedsHigherLimit
        logmsg(0,'** Spectrum exceeds frequency range. Artifacts at lower and upper frequency limits possible.');
      else
        logmsg(0,'** Spectrum exceeds frequency range. Artifacts at lower frequency limit possible.');
      end
    else
      if exceedsHigherLimit
        logmsg(0,'** Spectrum exceeds frequency range. Artifacts at upper frequency limit possible.');
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

  % Lorentzian broadening
  if fwhmL~=0
    if fwhmL>2*Exp.deltaX
      logmsg(1,'  convoluting with Lorentzian, FWHM %g MHz, derivative %d',fwhmL,HarmonicL);
      if min(size(spec))==1, fwhm = [fwhmL 0]; else, fwhm = [0 fwhmL]; end
      spec = convspec(spec,Exp.deltaX,fwhm,HarmonicL,0);
    else
      if HarmonicL==0
        % Skip convolution, since it has no effect with such a narrow delta-like Lorentzian.
      else
        error('Lorentzian linewidth is smaller than increment - cannot perform convolution.');
      end
    end
  end
  
  % Gaussian broadening
  if fwhmG~=0
    if fwhmG>2*Exp.deltaX
      logmsg(1,'  convoluting with Gaussian, FWHM %g MHz, derivative %d',fwhmG,HarmonicG);
      if min(size(spec))==1, fwhm = [fwhmG 0]; else, fwhm = [0 fwhmG]; end
      spec = convspec(spec,Exp.deltaX,fwhm,HarmonicG,1);
    else
      if HarmonicG==0
        % Skip convolution, since it has no effect with such a narrow delta-like Gaussian.
      else
        error('Gaussian linewidth is smaller than increment - cannot perform convolution.');
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
  logmsg(1,'  no linewidth given: returning stick spectrum');
  
  if Exp.Harmonic>0
    logmsg(1,'  harmonic %d: using differentiation',Exp.Harmonic);
    for h = 1:Exp.Harmonic
      dspec = diff(spec,[],2)/Exp.deltaX;
      spec = (dspec(:,[1 1:end]) + dspec(:,[1:end end]))/2;
    end
  else
    logmsg(1,'  harmonic 0: returning absorption spectrum');
  end

end

% Assign output.
%-----------------------------------------------------------------------
switch nargout
  case 0
    cla
    plot(xAxis,spec);
    axis tight
    xlabel('frequency (MHz)');
    ylabel('intensity (arb.u.)');
    if ~isnan(Exp.mwFreq)
      title(sprintf('ENDOR at %0.8g GHz and %0.8g mT',Exp.mwFreq,Exp.Field));
    else
      title(sprintf('ENDOR at %0.8g mT',Exp.Field));
    end
  case 1
    varargout = {spec};
  case 2
    varargout = {xAxis,spec};
  case 3
    info.Transitions = Transitions;
    varargout = {xAxis,spec,info};
end

% Report performance.
%-----------------------------------------------------------------------
hmsString = elapsedtime(StartTime,clock);
logmsg(1,['salt took ' hmsString]);

logmsg(1,'=end=salt=========%s=================\n',char(datetime));

clear global EasySpinLogLevel

end


function salt_plot(xAxis,spec,Exp)
cla
plot(xAxis,spec);
axis tight
xlabel('frequency (MHz)');
ylabel('intensity (arb.u.)');
if isfield(Exp,'mwFreq') && ~isnan(Exp.mwFreq)
  title(sprintf('ENDOR at %0.8g GHz and %0.8g mT',Exp.mwFreq,Exp.Field));
else
  title(sprintf('ENDOR at %0.8g mT',Exp.Field));
end
end
