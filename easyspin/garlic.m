% garlic    Simulates isotropic and fast-motion cw EPR spectra 
%
%   garlic(Sys,Exp)
%   garlic(Sys,Exp,Opt)
%   spec = ...
%   [B,spec] = ...
%
%   Computes the solution cw EPR spectrum of systems with
%   an unpaired electron and arbitrary numbers of nuclear spins.
%
%   Sys: spin system structure
%     g            isotropic g factor or 3 principal values of g
%     Nucs         string with comma-separated list of isotopes
%     n            vector of number of equivalent nuclei (default all 1)
%     A            vector of hyperfine couplings (MHz)
%     lw           vector with FWHM line widths (mT)
%                   1 element:  GaussianFWHM
%                   2 elements: [GaussianFWHM LorentzianFWHM]
%     lwpp         peak-to-peak line widths (mT), same format as Sys.lw
%
%     tcorr        correlation time for fast-motion linewidths (s)
%                  If omitted or zero, the isotropic spectrum is computed.
%     logtcorr     log10 of the correlation time for fast-motion linewidths (s)
%                  If logtcorr is given, tcorr is ignored.
%
%   Exp:  experimental parameter settings
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
%  Opt:  simulation parameters
%      Verbosity    log level (0 none, 1 normal, 2 very verbose)
%      Method       method used to calculate line positions:
%                     'exact' - Breit-Rabi solver
%                     'perturb1','perturb2','perturb3','perturb4','perturb5'
%      AccumMethod  method used to construct spectrum
%                     'binning' - stick spectrum with line shape convolution
%                     'template' - interpolative construction using line shape template
%                     'explicit' - explicit evaluation of line shape for each line
%      IsoCutoff    relative isotopologue abundance cutoff threshold
%                     between 0 and 1, default 1e-6
%
%   Output
%     B                magnetic field axis (mT)
%     spec             spectrum (arbitrary units)
%
%     If no output parameter is specified, the simulated spectrum
%     is plotted.

function varargout = garlic(Sys,Exp,Opt)

varargout = cell(1,nargout);
if nargin==0, help(mfilename); return; end

% Check expiry date
error(eschecker);

% Check Matlab version.
error(chkmlver);

switch nargin
  case 1, error('Experimental parameters (2nd input) are missing!');
  case 2, Opt = [];
  case 3
otherwise, error('Wrong number of input parameters!');
end

switch nargout
case {0,1,2,3}
otherwise, error('Wrong number of output parameters!');
end

if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 0; % Log level
end

global EasySpinLogLevel
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

if ~isfield(Opt,'IsoCutoff')
  Opt.IsoCutoff = 1e-6;
else
  if Opt.IsoCutoff<0 || Opt.IsoCutoff>1
    error('Options.IsoCutoff must be between 0 and 1.');
  end
end

if ~isfield(Opt,'Output'), Opt.Output = 'summed'; end
[Output,err] = parseoption(Opt,'Output',{'summed','separate'});
error(err);
summedOutput = Output==1;

if ~isfield(Sys,'singleiso') || ~Sys.singleiso
  
  if ~iscell(Sys), Sys = {Sys}; end
  
  nComponents = numel(Sys);
  if nComponents>1
    logmsg(1,'  %d component spin systems...');
  else
    logmsg(1,'  single spin system');
  end
  
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
      [xAxis,spec_,Transitions] = garlic(Sys_,Exp,Opt);
      
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
    case 1, varargout = {spec};
    case 2, varargout = {xAxis,spec};
    case 3, varargout = {xAxis,spec,Transitions};
  end
  return
end
%==================================================================

if EasySpinLogLevel>=1
  logmsg(1,['=begin=garlic=====' datestr(now) '=================']);
end

%-------------------------------------------------------------------------
% System structure check
%-------------------------------------------------------------------------

[Sys,err] = validatespinsys(Sys);
error(err);
if Sys.MO_present, error('salt does not support Sys.Ham* parameters!'); end
if any(Sys.L(:)), error('salt does not support Sys.L!'); end

if (Sys.nElectrons~=1) || (isfield(Sys,'L') && any(Sys.L))
  error('Only systems with one electron spin S=1/2 are supported.');
elseif Sys.S~=1/2
  error('Only systems with electron spin S=1/2 are supported.');
end
if isfield(Sys,'Exchange') && (Sys.Exchange~=0)
  error('garlic does not support Heisenberg exchange.');
end
if isfield(Sys,'nn') && any(Sys.nn(:)~=0)
  error('garlic does not support nuclear-nuclear couplings (Sys.nn).');
end

if isfield(Sys,'logtcorr'), Sys.tcorr = 10.^Sys.logtcorr; end
if ~isfield(Sys,'tcorr'), Sys.tcorr = 0; end

FastMotionRegime = ~isempty(Sys.tcorr) && (Sys.tcorr~=0);
if FastMotionRegime
  if ~isreal(Sys.tcorr)
    error('Problem with System.tcorr: must be a number!');
  end
  if numel(Sys.tcorr)~=1
    error('Problem with System.tcorr: must be a single number!');
  end
  if Sys.tcorr<=0
    error('Problem with System.tcorr: must be positive!');
  end
  if Sys.tcorr<1e-13
    error('Correlation time too small (%f seconds).',Sys.tcorr);
  end
  if isfield(Sys,'n')
    if any(Sys.n>1)
      error('Cannot treat equivalent nuclei in fast-motion regime!\n Please rewrite spin system!');
    end
  end
  if isfield(Sys,'Potential')
    error('garlic cannot simulate fast-motion spectra with an orienting potential (Sys.Potential).');
  end
end

ConvolutionBroadening = any(Sys.lw>0) || FastMotionRegime;

%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% Checking experiment structure
%-------------------------------------------------------------------------
DefaultExp.Harmonic = [];
DefaultExp.nPoints = 1024;
DefaultExp.ModAmp = 0;
DefaultExp.Mode = 'perpendicular';
DefaultExp.mwPhase = 0;
DefaultExp.Temperature = NaN; % don't compute thermal equilibrium polarizations

Exp = adddefaults(Exp,DefaultExp);

% Microwave frequency
if ~isfield(Exp,'mwFreq') || isempty(Exp.mwFreq)
  if ~isfield(Exp,'Field')
    error('Please supply either the microwave frequency in Exp.mwFreq (for field sweeps) or the magnetic field in Exp.Field (for frequency sweeps).');
  end
  FieldSweep = false;
else
  if isfield(Exp,'Field')
    if ~isempty(Exp.Field)
      error('Give either Exp.mwFreq (for a field sweep) or Exp.Frequency (for a frequency sweep), but not both.');
    end
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

% Sweep range (magnetic field or frequency)
SweepAutoRange = false;
if FieldSweep
  if isfield(Exp,'CenterSweep')
    if isfield(Exp,'Range')
      logmsg(0,'Using Exp.CenterSweep and ignoring Exp.Range.');
    end
  else
    if isfield(Exp,'Range')
      if (Exp.Range(1)>=Exp.Range(2)) || any(Exp.Range<0)
        error('Invalid sweep range!');
      end
      Exp.CenterSweep = [mean(Exp.Range) diff(Exp.Range)];
    else
      logmsg(1,'  automatic determination of sweep range');
      SweepAutoRange = true;
    end
  end
  if ~SweepAutoRange
    Exp.Range = Exp.CenterSweep(1) + [-1 1]/2*Exp.CenterSweep(2);
  end
else
  if isfield(Exp,'mwCenterSweep')
    if isfield(Exp,'mwRange')
      logmsg(0,'Using Exp.mwCenterSweep and ignoring Exp.mwRange.');
    end
  else
    if isfield(Exp,'mwRange')
      if (Exp.mwRange(1)>=Exp.mwRange(2)) || any(Exp.mwRange<0)
        error('Invalid sweep range!');
      end
      Exp.mwCenterSweep = [mean(Exp.mwRange) diff(Exp.mwRange)];
    else
      logmsg(1,'  automatic determination of sweep range');
      SweepAutoRange = true;
    end
  end
  if ~SweepAutoRange
    Exp.mwRange = Exp.mwCenterSweep(1) + [-1 1]/2*Exp.mwCenterSweep(2);
  end
end

% Number of points
if any(~isreal(Exp.nPoints)) || numel(Exp.nPoints)>1 || (Exp.nPoints<2)
  error('Problem with Exp.nPoints. Needs to be a number not smaller than 2.')
end

% Detection harmonic
autoHarmonic = ~isfield(Exp,'Harmonic') || isempty(Exp.Harmonic) || isnan(Exp.Harmonic);
noBroadening = ~FastMotionRegime && all(Sys.lw==0) && all(Sys.lwpp==0);
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

% Resonator mode
if strcmp(Exp.Mode,'perpendicular')
  ParallelMode = 0;
elseif strcmp(Exp.Mode,'parallel')
  ParallelMode = 1;
else
  error('Exp.Mode must be either ''perpendicular'' or ''parallel''.');
end
logmsg(1,'  harmonic %d, %s mode',Exp.Harmonic,Exp.Mode);

% Temperature
if ~isnan(Exp.Temperature)
  if isinf(Exp.Temperature)
    error('If given, Exp.Temperature must be a finite value.')
  end
  if numel(Exp.Temperature)~=1
    error('garlic does not support spin polarization in Exp.Temperature.');
  end
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


% Complain if fields only valid in pepper() are given
if isfield(Exp,'Orientations')
  warning('Exp.Orientations is not used by garlic.');
end
if isfield(Exp,'Ordering')
  warning('Exp.Ordering is not used by garlic.');
end
if isfield(Exp,'CrystalSymmetry')
  warning('Exp.CrystalSymmetry is not used by garlic.');
end
if isfield(Exp,'CrystalOrientation')
  warning('Exp.CrystalOrientation is not used by garlic.');
end


%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% Simulation options
%-------------------------------------------------------------------------

% Stretch factor for automatically detected field range
if ~isfield(Opt,'Stretch')
  Opt.Stretch = 0.25;
end

% Maximum number of Newton-Raphson or fixed-point iteration steps
if ~isfield(Opt,'MaxIterations')
  Opt.MaxIterations = 15;
end

% Accuracy of resonance field computation:
% If relative iteration-to-iteration change of resonance field
% falls below this value, the computation is stopped.
if ~isfield(Opt,'Accuracy')
  Opt.Accuracy = 1e-12;
  % Note: relative accuracy of Planck constant is 5e-8
end

%Options.AcceptLimit = 0.5; % rejection threshold for large hyperfine couplings

% Method for computation of the resonance fields
if isfield(Opt,'PerturbOrder')
  error('Options.PerturbOrder is obsolete. Use Options.Method = ''perturb2'' etc instead.');
end
if ~isfield(Opt,'Method'), Opt.Method = 'exact'; end
switch Opt.Method
  case 'exact',    PerturbOrder = 0; % exact calculation, Breit/Rabi
  case 'perturb',  PerturbOrder = 5;
  case 'perturb1', PerturbOrder = 1;
  case 'perturb2', PerturbOrder = 2;
  case 'perturb3', PerturbOrder = 3;
  case 'perturb4', PerturbOrder = 4;
  case 'perturb5', PerturbOrder = 5;
  otherwise
    error('Unknown method ''%s''.',Opt.Method);
end

% Method for spectrum construction
if ~isfield(Opt,'AccumMethod') || isempty(Opt.AccumMethod)
  if FastMotionRegime
    Opt.AccumMethod = 'explicit';
  else
    Opt.AccumMethod = 'binning';
  end
end

%-------------------------------------------------------------------------

if Sys.fullg
  giso = mean(eig(Sys.g(1:3,:)));
else
  giso = mean(Sys.g(1,:));
end
if FieldSweep
  mwFreq = Exp.mwFreq*1e9; % GHz -> Hz
  CentralResonance = (mwFreq*planck)/(bmagn*giso)*1e3; % mT, absence of nuclei
else
  CentralResonance = bmagn*giso*Exp.Field*1e-3/planck; % Hz
end

if FastMotionRegime
  if FieldSweep
    FastMotionLw = fastmotion(Sys,CentralResonance,Sys.tcorr,'field'); % mT
  else
    FastMotionLw = fastmotion(Sys,Exp.Field,Sys.tcorr,'freq')/1e3; % MHz -> GHz
  end
  if all(FastMotionLw==0)
    error('Linewidths resulting from fast-motion lindwidth parameters must be positive! Did you supply isotropic values only?');
  elseif any(FastMotionLw<=0)
    error('Linewidths resulting from fast-motion lindwidth parameters must be positive!');
  end
end

if Sys.nNuclei>0
  I_all = Sys.I;
  gn = Sys.gn;
  if ~Sys.fullA
    a_all = mean(Sys.A,2).'*1e6; % isotropic A values, MHz -> Hz
  else
    for iNuc=1:Sys.nNuclei
      Adiag_ = eig(Sys.A((iNuc-1)*3+(1:3),:));
      a_all(iNuc) = mean(Adiag_)*1e6; % isotropic A values, MHz -> Hz
    end
  end
  n_all = Sys.n;
  if FieldSweep
    % Make sure mw frequency is larger than zero field splitting.
    zfSplitting = sum(abs(a_all).*(I_all+0.5)); % approximate
    if mwFreq<zfSplitting
      fprintf(['Microwave frequency (%f GHz) is smaller than\nlargest splitting ' ...
        'of levels at zero field (estimate %f GHz).\nSpectrum might be inaccurate.\n' ...
        'Consider using the function pepper() instead.\n'],mwFreq/1e9,zfSplitting/1e9);
    end
  end
  if any(I_all==0)
    error('Nuclei with spin 0 present. Cannot compute spectrum.');
  end
  if any(a_all==0)
    error('Nuclei with hyperfine coupling 0 present. Cannot compute spectrum.');
  end
  a_all = a_all*planck;   % Hz -> Joule
end

logmsg(1,'Computing resonance shifts...');

Shifts = cell(1,Sys.nNuclei);
Amplitudes = cell(1,Sys.nNuclei);

for iNucGrp = 1:Sys.nNuclei
  
  if n_all(iNucGrp)==0
    Shifts{iNucGrp} = 0;
    Amplitudes{iNucGrp} = 1;
    continue
  end

  [I_this,nn] = equivcouple(I_all(iNucGrp),n_all(iNucGrp));
  nLines = (2*I_all(iNucGrp)+1)^n_all(iNucGrp);
  nFSpins = numel(nn);
  aiso = a_all(iNucGrp); % Joule
  aiso = abs(aiso);
  
  Positions = [];
  Intensities = [];

  if FieldSweep
    
    % Field sweep
    %-------------------------------------------------------------------
    if PerturbOrder==0
      logmsg(1,'  Breit-Rabi solver, accuracy %g',Opt.Accuracy);
      maxIterationsDone = -1;
      
      % Fixed-point iteration, based on the Breit-Rabi formula
      % J. A. Weil, J. Magn. Reson. 4, 394-399 (1971)
      gammae = giso*bmagn;
      h = planck;
      for iF = 1:nFSpins
        if nn(iF)==0, continue; end
        I = I_this(iF);
        mI = -I:I;
        B_ = 0;
        RelChangeB = inf;
        gamman = 0; % to obtain start field, neglect NZI
        for iter = 1:Opt.MaxIterations
          nu = mwFreq + gamman/h*B_;
          q = 1 - (aiso./(2*h*nu)).^2;
          Bnew = aiso./(gammae+gamman)./q.* ...
            (-mI+sign(aiso)*sqrt(mI.^2+q.*((h*nu/aiso).^2 - (I+1/2)^2)));
          if ~isreal(Bnew)
            error('Hyperfine coupling too large compared to microwave frequency. Cannot compute resonance field positions. Use pepper() instead.');
          end
          if B_~=0, RelChangeB = 1 - Bnew./B_; end
          B_ = Bnew;
          gamman = gn(iNucGrp)*nmagn; % re-include NZI
          if abs(RelChangeB)<Opt.Accuracy, break; end
        end
        if iter>maxIterationsDone, maxIterationsDone = iter; end
        if iter==Opt.MaxIterations
          error(sprintf('Breit-Rabi solver didn''t converge after %d iterations!',iter));
        end
        Positions = [Positions B_];
        Intensities = [Intensities nn(iF)/nLines*ones(size(B_))];
      end
      logmsg(1,'  maximum %d iterations done',maxIterationsDone);
      
    else
      logmsg(1,'  perturbation expansion, order %d, accuracy %g',PerturbOrder,Opt.Accuracy);
      
      % Based on \Delta E - h nu = 0 with a 5th-order Taylor expansion in aiso of
      % the Breit-Rabi expression for \Delta E. The resulting polynomial in aiso with terms
      % (0,aiso,aiso^2,aiso^3,aiso^4,...) is also a polynomial in B with the same terms
      % (B,0,1/B,1/B^2,1/B^3,...)
      maxIterationsDone = -1;
      pre = aiso .* (aiso/(giso*bmagn+gn(iNucGrp)*nmagn)/2).^(1:4);
      for iF = 1:nFSpins
        if nn(iF)==0, continue; end
        I = I_this(iF);
        for mI = -I:I
          c(6) = pre(4) * mI * (1+6*(I-1)*I*(I+1)*(I+2)-10*(-3+2*I*(I+1))*mI^2+14*mI^4);
          c(5) = pre(3) *      (I-2*I^3-I^4+6*(I^2+I-1)*mI^2-5*mI^4);
          c(4) = pre(2) * mI * (1-2*I*(I+1)+2*mI^2);
          c(3) = pre(1) *      (I*(I+1)-mI^2);
          c(2) = aiso*mI - mwFreq*planck;
          c(1) = giso*bmagn;
          % Newton-Raphson
          % Usually between 2 and 4 steps are necessary, never more than 7
          n = PerturbOrder+1;
          c = c(1:n);
          dc = (n-1:-1:1).*c(1:n-1);
          Bb = -c(2)/c(1); % first-order as start guess
          dB = Bb;
          for iter = 1:Opt.MaxIterations
            if abs(dB/Bb)<Opt.Accuracy, break; end
            dB = polyval(c,Bb)/polyval(dc,Bb);
            Bb = Bb - dB;
          end
          if iter>maxIterationsDone, maxIterationsDone = iter; end
          if iter==Opt.MaxIterations
            error('Newton-Raphson has convergence problem!');
          end
          Positions(end+1) = Bb; % T
          Intensities(end+1) = nn(iF)/nLines;
        end
      end
      logmsg(1,'  maximum %d iterations done',maxIterationsDone);
      
    end
    Positions = Positions*1e3; % T -> mT

  else
    
    % Frequency sweep
    %-------------------------------------------------------------------
    if PerturbOrder==0
      logmsg(1,'  Breit-Rabi formula');
      
      % Breit-Rabi formula, J. A. Weil, J. Magn. Reson. 4, 394-399 (1971), [1]
      B_ = Exp.Field*1e-3; % T
      gebB = giso*bmagn*B_;
      gnbB = gn(iNucGrp)*nmagn*B_;
      for iF = 1:nFSpins
        if nn(iF)==0, continue; end
        I = I_this(iF);
        mI = -I:I;
        alpha = (gebB+gnbB)/aiso/(I+1/2);
        E1 = -aiso/4-gnbB*(mI+1/2)+(I+1/2)*aiso/2*sqrt(1+2*(mI+1/2)/(I+1/2)*alpha+alpha^2);
        E2 = -aiso/4-gnbB*(mI-1/2)-(I+1/2)*aiso/2*sqrt(1+2*(mI-1/2)/(I+1/2)*alpha+alpha^2);
        nu_ = (E1-E2)/planck; % J -> Hz
        Positions = [Positions nu_]; % Hz
        Intensities = [Intensities nn(iF)/nLines*ones(size(nu_))];
      end
      
    else
      logmsg(1,'  perturbation theory, order %d',PerturbOrder);
      
      % Taylor expansion in aiso around 0 of the Breit-Rabi expression for \Delta E.
      B_ = Exp.Field*1e-3;
      pre = aiso .* (aiso/(giso*bmagn*B_+gn(iNucGrp)*nmagn*B_)/2).^(1:4);
      for iF = 1:nFSpins
        if nn(iF)==0, continue; end
        I = I_this(iF);
        for mI = -I:I
          c(1) = giso*bmagn*B_; % all c(:) in Joule
          c(2) = aiso*mI;
          c(3) = pre(1) *      (I*(I+1)-mI^2);
          c(4) = pre(2) * mI * (1-2*I*(I+1)+2*mI^2);
          c(5) = pre(3) *      (I-2*I^3-I^4+6*(I^2+I-1)*mI^2-5*mI^4);
          c(6) = pre(4) * mI * (1+6*(I-1)*I*(I+1)*(I+2)-10*(-3+2*I*(I+1))*mI^2+14*mI^4);
          dE = sum(c(1:PerturbOrder+1));
          nu_ = dE/planck; % Joule -> Hz
          Positions(end+1) = nu_; % Hz
          Intensities(end+1) = nn(iF)/nLines;
        end
      end
      
    end

  end
  Shifts{iNucGrp} = Positions - CentralResonance; % field sweep: mT; freq sweep: Hz
  Amplitudes{iNucGrp} = Intensities;
  
  logmsg(1,'  spin group %d: %d F spins, %d lines',iNucGrp,nFSpins,numel(Shifts{iNucGrp}));
end

% Statistics: minimum, maximum, number of lines
%--------------------------------------------------------------
nPeaks = 1;
posmax = CentralResonance;
posmin = CentralResonance;
for iShift = 1:numel(Shifts)
  nPeaks = nPeaks * numel(Shifts{iShift});
  posmax = posmax + max(Shifts{iShift});
  posmin = posmin + min(Shifts{iShift});
end
logmsg(1,'  total %d lines',nPeaks);


if FastMotionRegime
  % Add Lorentzian
  LorentzianLw = FastMotionLw + Sys.lw(2);
  Sys.lw(2) = 0;
  maxLw = max(max(FastMotionLw),Sys.lw(1));
else
  LorentzianLw = Sys.lw(2);
  maxLw = max(Sys.lw);
end

% Autoranging
%--------------------------------------------------------------
if SweepAutoRange
  if FieldSweep
    posrange = (posmax-posmin)*Opt.Stretch; % mT
    minrange = 1; % mT
    Exp.Range = [posmin,posmax] + [-1 1]*max([5*maxLw,posrange,minrange]);
    Exp.Range(1) = max(Exp.Range(1),0);
    logmsg(1,'  automatic field range from %g mT to %g mT',Exp.Range(1),Exp.Range(2));
  else
    posrange = (posmax-posmin)*Opt.Stretch; % Hz
    minrange = 1e6; % Hz
    Exp.mwRange = [posmin,posmax] + [-1 1]*max([5*maxLw/1e3,posrange,minrange]);
    Exp.mwRange(1) = max(Exp.mwRange(1),0);
    Exp.mwRange = Exp.mwRange/1e9; % Hz -> GHz
    logmsg(1,'  automatic frequency range from %g GHz to %g GHz',Exp.mwRange(1),Exp.mwRange(2));
  end
end

if FieldSweep
  logmsg(1,'  spectral spread %g mT\n  g=%g resonance at %g mT',(posmax-posmin),giso,CentralResonance);
else
  logmsg(1,'  spectral spread %g MHz\n  g=%g resonance at %g GHz',(posmax-posmin)/1e6,giso,CentralResonance/1e9);
end

% Combining shifts and intensities
%--------------------------------------------------------------
logmsg(1,'Combining line shifts and line multiplicities...');
if nPeaks>1
  Positions = allcombinations(Shifts{:},'+') + CentralResonance;
  Intensities = allcombinations(Amplitudes{:},'*');
else
  Positions = CentralResonance;
  Intensities = 1;
end
if ~FieldSweep
  Positions = Positions/1e9; % Hz -> GHz
end

% Line intensities
%--------------------------------------------------------------
logmsg(1,'Computing overall line intensities...');
if ParallelMode
  % Parallel mode: no intensities
  Intensities = zeros(size(Intensities));
else

  % Transition rate
  g1mean2 = giso^2;
  TransitionRate = (8*pi^2)*g1mean2*(bmagn/planck/1e9/2)^2;
  Intensities = Intensities*TransitionRate;
  
  % 1/g factor (mT/MHz)
  if FieldSweep
    dBdE = planck/(giso*bmagn)*1e9;
    Intensities = Intensities*dBdE;
  end
  
  % Temperature: thermal equilibrium polarization
  if isfinite(Exp.Temperature)
    if FieldSweep
      DeltaE = planck*Exp.mwFreq*1e9;  % Joule
    else
      DeltaE = bmagn*giso*Exp.Field*1e-3;  % Joule
    end
    e = exp(-DeltaE/boltzm/Exp.Temperature);
    Population = [1 e]/(1+e);
    Polarization = Population(1) - Population(2);
    Intensities = Intensities*Polarization;
  end
  
end


%===================================================================
% Spectrum construction
%===================================================================

% (1) Initial spectrum construction
%-------------------------------------------------------------------
if FieldSweep
  SweepRange = Exp.Range;
else
  SweepRange = Exp.mwRange;
end
xAxis = linspace(SweepRange(1),SweepRange(2),Exp.nPoints);

switch Opt.AccumMethod
  case 'template'
    % Accumulate spectrum by interpolating from a pre-computed template lineshape
    
    logmsg(1,'Constructing spectrum using Lorentzian lineshape template...');
  
    if LorentzianLw==0
      error('Cannot use templated linshape accumulation with zero linewidth.');
    end
    if Exp.mwPhase~=0
      error('Cannot use templated linshape accumulation with non-zero Exp.mwPhase.');
    end
    
    dxFine = min(xAxis(2)-xAxis(1),min(LorentzianLw)/5);
    nPointsFine = round((SweepRange(2)-SweepRange(1))/dxFine+1);
    xAxisFine = linspace(SweepRange(1),SweepRange(2),nPointsFine);
    dxFine = xAxisFine(2) - xAxisFine(1);
    
    xT = 1e5;
    wT = xT/20; % 0.0025 at borders for Harmonic = -1
    Template = lorentzian(0:2*xT-1,xT,wT,Exp.Harmonic-1);
    if numel(LorentzianLw)==1
      LorentzianLw = LorentzianLw*ones(size(Positions));
    end
    spec = lisum1i(Template,xT,wT,Positions,Intensities,LorentzianLw,xAxisFine);
    Exp.ConvHarmonic = 0;
    LorentzianLw = 0;

  case 'explicit'
    % Accumulate spectrum by explicit evaluation of lineshape function
    
    logmsg(1,'Constructing spectrum with explicit Lorentzian lineshapes...');
    if LorentzianLw==0
      error('Cannot use axplicit accumulation with zero Lorentzian linewidth.');
    end
    
    dxFine = min(xAxis(2)-xAxis(1),min(LorentzianLw)/5);
    nPointsFine = round((SweepRange(2)-SweepRange(1))/dxFine+1);
    xAxisFine = linspace(SweepRange(1),SweepRange(2),nPointsFine);
    dxFine = xAxisFine(2) - xAxisFine(1);
    
    spec = 0;
    if numel(LorentzianLw)==1
      for iLine = 1:numel(Positions)
        spec = spec + Intensities(iLine)*lorentzian(xAxisFine,Positions(iLine),LorentzianLw,Exp.ConvHarmonic,Exp.mwPhase);
      end
    else
      for iLine = 1:numel(Positions)
        spec = spec + Intensities(iLine)*lorentzian(xAxisFine,Positions(iLine),LorentzianLw(iLine),Exp.ConvHarmonic,Exp.mwPhase);
      end
    end
    Exp.mwPhase = 0;
    Exp.ConvHarmonic = 0;
    LorentzianLw = 0;
    
  case 'binning'
    % Accumulate spectrum by binning of delta peaks
    
    logmsg(1,'Constructing stick spectrum using binning...');
    
    if FastMotionRegime
      error('Cannot use delta binning (Options.AccumMethod=''binning'' in the fast-motion regime.');
    end
    
    expandFactor = 1;
    nPointsFine = (Exp.nPoints-1)*expandFactor + 1;
    xAxisFine = linspace(SweepRange(1),SweepRange(2),nPointsFine);
    dxFine = xAxisFine(2) - xAxisFine(1);
    
    spec = constructstickspectrum(Positions,Intensities,SweepRange,nPointsFine);
    spec = spec/dxFine;
    
  otherwise
    error('\n  Unknown method ''%s'' in Options.AccumMethod.\n',Opt.AccumMethod);

end

% (2) Convolutional broadening
%-------------------------------------------------------------------
if ConvolutionBroadening
  
  fwhmL = LorentzianLw;
  fwhmG = Sys.lw(1);
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
    unitstr = 'MHz';
    fwhmL = fwhmL/1e3; % MHz -> GHz
    fwhmG = fwhmG/1e3; % MHz -> GHz
  end
  
  % Convolution with Lorentzian
  if fwhmL~=0
    if fwhmL>2*dxFine
      logmsg(1,'Convoluting with Lorentzian (FWHM %g %s, derivative %d)...',fwhmL,unitstr,HarmonicL);
      spec = convspec(spec,dxFine,fwhmL,HarmonicL,0,mwPhaseL);
    else
      if HarmonicL==0
        % Skip convolution, since it has no effect with such a narrow delta-like Lorentzian.
      else
        error('Lorentzian linewidth is smaller than increment - cannot perform convolution.');
      end
    end
  end
  
  % Convolution with Gaussian
  if fwhmG~=0
    if fwhmG>2*dxFine
      logmsg(1,'Convoluting with Gaussian (FWHM %g %s, derivative %d)...',fwhmG,unitstr,HarmonicG);
      spec = convspec(spec,dxFine,fwhmG,HarmonicG,1,mwPhaseG);
    else
      if HarmonicG==0
        % Skip convolution, since it has no effect with such a narrow delta-like Gaussian.
      else
        error('Gaussian linewidth is smaller than increment - cannot perform convolution.');
      end
    end
  end
  
end

if numel(xAxisFine)~=numel(xAxis)
  logmsg(1,'Re-interpolation (%d -> %d points)...',numel(xAxisFine),numel(xAxis));
  spec = interp1(xAxisFine,spec,xAxis);
end

% (3) Field modulation
%-------------------------------------------------------------------
if FieldSweep
  if Exp.ModAmp>0
    logmsg(1,'Applying field modulation (%g mT amplitude)...',Exp.ModAmp);
    spec = fieldmod(xAxis,spec,Exp.ModAmp,Exp.ModHarmonic);
  else
    % derivatives already included in the convolution
  end
else
  % frequency sweeps: modulation not available
end

%===================================================================

switch nargout
  case 1, varargout = {spec};
  case 2, varargout = {xAxis,spec};
  case 3, varargout = {xAxis,spec,Positions};
end
if EasySpinLogLevel>=1
  logmsg(1,'=end=garlic=======%s=================\n',datestr(now));
end
clear global EasySpinLogLevel

return


%===================================================================
%===================================================================
%===================================================================
function Spectrum = constructstickspectrum(Positions,Amplitudes,Range,nPoints)

% Convert positions to indices into spectral vector (1...nPoints)
idxPositions = (Positions-Range(1))/diff(Range) * (nPoints-1);
idxPositions = 1 + round(idxPositions);

% Identify in-range lines
inRange = (idxPositions>=1) & (idxPositions<=nPoints);
if any(idxPositions<1)  
  logmsg(0,'** Spectrum exceeds sweep range. Artifacts at lower limit possible.');
end
if any(idxPositions>nPoints)
  logmsg(0,'** Spectrum exceeds sweep range. Artifacts at upper limit possible.');
end

% Bin in-range lines into spectrum
Spectrum = full(sparse(1,idxPositions(inRange),Amplitudes(inRange),1,nPoints));

return
