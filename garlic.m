% garlic    Simulates isotropic and fast-motion cw EPR spectra 
%
%   garlic(Sys,Exp)
%   spec = garlic(Sys,Exp)
%   [B,spec] = garlic(Sys,Exp)
%
%   Computes the solution cw EPR spectrum of systems with
%   an unpaired electron and arbitrary numbers of nuclear spins.
%
%   Input:
%     Sys.g            isotropic g factor or 3 principal values of g
%     Sys.Nucs         string with comma-separated list of isotopes
%     Sys.n            vector of number of equivalent nuclei (default all 1)
%     Sys.A            vector of hyperfine couplings [MHz]
%     Sys.lw           vector with FWHM line widths [mT]
%                       1 element:  GaussianFWHM
%                       2 elements: [GaussianFWHM LorentzianFWHM]
%     Sys.lwpp         peak-to-peak line widths [mT], same format as Sys.lw
%
%     Sys.tcorr        correlation time for fast-motion linewidths [s]
%                      If omitted or zero, the isotropic spectrum is computed.
%     Sys.logtcorr     log10 of the correlation time for fast-motion linewidths [s]
%                      If logtcorr is given, tcorr is ignored.
%
%     Exp.mwFreq       spectrometer frequency [GHz]
%     Exp.CenterSweep  [centerField sweepRange] in mT
%     Exp.Range        [minField maxField] in mT
%          Exp.Range is only used if Exp.CenterSweep is not given.
%          If both Exp.CenterField and Exp.Range are omitted, the
%          magnetic field range is determined automatically.
%
%     Exp.Harmonic     detection harmonic (0, 1, or 2), default 1
%     Exp.nPoints      number of points (default 1024)
%
%   Output
%     B                magnetic field axis [mT]
%     spec             spectrum [arbitrary units]
%
%     If no output parameter is specified, the simulated spectrum
%     is plotted.
%
%   Example isotropic spectrum
%
%     Sys = struct('g',2,'Nucs','1H,14N','A',[30,40],'n',[3,4]);
%     Sys.lwpp = [0,0.1];
%     Exp = struct('mwFreq',9.7,'nPoints',10000);
%     garlic(Sys,Exp);
%
%   Example fast-motion spectrum
%
%     A = [16, 16, 86];
%     Sys = struct('g',[2.0088 2.0061 2.0027],'Nucs','14N','A',A);
%     Sys.tcorr = 1e-9;
%     Exp = struct('mwFreq',9.5);
%     garlic(Sys,Exp);

function varargout = garlic(Sys,Exp,Options)

varargout = cell(1,nargout);
if (nargin==0), help(mfilename); return; end

switch (nargin)
  case 1, error('Experimental parameters (2nd input) are missing!');
  case 2, Options = [];
  case 3,
otherwise, error('Wrong number of input parameters!');
end

switch (nargout)
case {0,1,2,3},
otherwise, error('Wrong number of output parameters!');
end

if ~isfield(Options,'Verbosity')
  Options.Verbosity = 0; % Log level
end

global EasySpinLogLevel;
EasySpinLogLevel = Options.Verbosity;

% Check Matlab version.
error(chkmlver);

% --------License ------------------------------------------------
LicErr = 'Could not determine license.';
Link = 'epr@eth'; eschecker; error(LicErr); clear Link LicErr
% --------License ------------------------------------------------


%==================================================================
% Loop over components and isotopologues
%==================================================================
FrequencySweep = ~isfield(Exp,'mwFreq') & isfield(Exp,'Field');

SweepAutoRange = (~isfield(Exp,'Range') || isempty(Exp.Range)) && ...
                 (~isfield(Exp,'CenterSweep') || isempty(Exp.CenterSweep));
if ~isfield(Options,'IsoCutoff'), Options.IsoCutoff = 1e-6; end

if ~isfield(Sys,'singleiso') || (Sys.singleiso==0)
  
  [SysList,weight] = expandcomponents(Sys,Options.IsoCutoff);
    
  if (numel(SysList)>1) && SweepAutoRange
    error('Multiple components: Please specify magnetic field range manually using Exp.Range or Exp.CenterSweep.');
  end
  
  spec = 0;
  for iComponent = 1:numel(SysList)
    [xAxis,spec_,Bres] = garlic(SysList{iComponent},Exp,Options);
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
    case 3, varargout = {xAxis,spec,Bres};
  end
  return
end
%==================================================================

if (EasySpinLogLevel>=1)
  logmsg(1,['=begin=garlic=====' datestr(now) '=================']);
end

%-------------------------------------------------------------------------
% System structure check
%-------------------------------------------------------------------------

[Sys,err] = validatespinsys(Sys);
error(err);


if (Sys.nElectrons~=1)
  error('Only systems with one electron spin S=1/2 are supported.');
elseif (Sys.S~=1/2)
  error('Only systems with electron spin S=1/2 are supported.');
end
if isfield(Sys,'Exchange') && (Sys.Exchange~=0)
  error('garlic does not support Heisenberg exchange.');
end

if isfield(Sys,'logtcorr'), Sys.tcorr = 10.^Sys.logtcorr; end
if ~isfield(Sys,'tcorr'), Sys.tcorr = 0; end

FastMotionRegime = (Sys.tcorr~=0);
if FastMotionRegime
  if (~isreal(Sys.tcorr))
    error('Problem with System.tcorr: must be a number!');
  end
  if (numel(Sys.tcorr)~=1)
    error('Problem with System.tcorr: must be a single number!');
  end
  if (Sys.tcorr<=0)
    error('Problem with System.tcorr: must be positive!');
  end
  if (Sys.tcorr<1e-13)
    error('Correlation time too small (%f seconds).',Sys.tcorr);
  end
  if isfield(Sys,'n');
    if any(Sys.n>1)
      error('Cannot treat equivalent nuclei in fast-motion regime!\n Please rewrite spin system!');
    end
  end
end

ConvolutionBroadening = any(Sys.lw>0);
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

% Number of points
if any(~isreal(Exp.nPoints)) || numel(Exp.nPoints)>1 || (Exp.nPoints<2)
  error('Problem with Exp.nPoints. Needs to be a number not smaller than 2.')
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
noBroadening = all(Sys.lw==0) && all(Sys.lwpp==0);
if (Exp.Harmonic>0) && ~FastMotionRegime && noBroadening
  error(['No linewidth/broadening given. Cannot compute spectrum with Exp.Harmonic=%d.\n'...
    'Please specify a line broadening).'],Exp.Harmonic);
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
    Exp.Harmonic = Exp.Harmonic - 1;
  else
    error('Exp.ModAmp cannot be used with frequency sweeps.');
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


%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% Simulation options
%-------------------------------------------------------------------------

% Stretch factor for automatically detected field range
if ~isfield(Options,'Stretch');
  Options.Stretch = 0.25;
end

% Maximum number of Newton-Raphson or fixed-point iteration steps
if ~isfield(Options,'MaxIterations')
  Options.MaxIterations = 15;
end

% Accuracy of resonance field computation:
% If relative iteration-to-iteration change of resonance field
% falls below this value, the computation is stopped.
if ~isfield(Options,'Accuracy')
  Options.Accuracy = 1e-12;
  % Note: relative accuracy of Planck constant is 5e-8
end

%Options.AcceptLimit = 0.5; % rejection threshold for large hyperfine couplings

% Method for computation of the resonance fields
if isfield(Options,'PerturbOrder')
  error('Options.PerturbOrder is obsolete. Use Opt.Method = ''perturb2'' etc instead.');
end
if ~isfield(Options,'Method'), Options.Method = 'exact'; end
switch Options.Method
  case 'exact',    PerturbOrder = 0; % exact calculation, Breit/Rabi
  case 'perturb',  PerturbOrder = 5;
  case 'perturb1', PerturbOrder = 1;
  case 'perturb2', PerturbOrder = 2;
  case 'perturb3', PerturbOrder = 3;
  case 'perturb4', PerturbOrder = 4;
  case 'perturb5', PerturbOrder = 5;
  otherwise
    error('Unknown method ''%s''.',Options.Method);
end

%-------------------------------------------------------------------------

if Sys.fullg
  giso = mean(mean(Sys.g(1:3,:)));
else
  giso = mean(Sys.g(1,:));
end
if FieldSweep
  mwFreq = Exp.mwFreq*1e9; % GHz -> Hz
  CentralResonance = (mwFreq*planck)/(bmagn*giso)*1e3; % mT, absence of nuclei
else
  CentralResonance = bmagn*giso*Exp.Field*1e-3/planck; % Hz
end

if (FastMotionRegime)
  if FieldSweep
    FastMotionLw = fastmotion(Sys,CentralResonance,Sys.tcorr);
  else
    FastMotionLw = fastmotion(Sys,Exp.Field,Sys.tcorr);
  end
  if all(FastMotionLw==0)
    error('Linewidths resulting from fast-motion lindwidth parameters must be positive! Did you supply isotropic values only?');
  elseif any(FastMotionLw<=0)
    error('Linewidths resulting from fast-motion lindwidth parameters must be positive!');
  end
end

if (Sys.nNuclei>0)
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
  a_all = a_all.*Sys.Ascale;
  n_all = Sys.n;
  if FieldSweep
    % Make sure zero field splitting is smaller than mw frequency.
    if sum(abs(a_all).*I_all)*2>mwFreq
      disp('Microwave frequency is smaller than largest splitting of levels at zero field. Spectrum might be inaccurate. Consider using the function pepper instead.');
    end
  end
  if any(I_all==0)
    error('Nuclei with spin 0 present.');
  end
  if any(a_all==0)
    error('Nuclei with coupling 0 present.');
  end
end

logmsg(1,'Computing resonance shifts...');
p = 1;
Shifts = {};
Amp = {};

for iNucGrp = 1:Sys.nNuclei
  if (n_all(iNucGrp)==0), continue; end

  [I_this,nn] = equivcouple(I_all(iNucGrp),n_all(iNucGrp));
  nLines = (2*I_all(iNucGrp)+1)^n_all(iNucGrp);
  nFSpins = length(nn);
  aiso = a_all(iNucGrp)*planck; % J
  Positions = [];
  Intensities = [];

  if FieldSweep
    
    % Field sweep
    %-------------------------------------------------------------------
    if (PerturbOrder==0)
      logmsg(1,'  Breit-Rabi solver, accuracy %g',Options.Accuracy);
      maxIterationsDone = -1;
      
      % Fixed-point iteration, based on the Breit-Rabi formula
      % J. A. Weil, J. Magn. Reson. 4, 394-399 (1971)
      gammae = giso*bmagn;
      h = planck;
      for iF = 1:nFSpins
        if (nn(iF)==0), continue; end
        I = I_this(iF);
        mI = -I:I;
        B_ = 0;
        RelChangeB = inf;
        gamman = 0; % to obtain start field, neglect NZI
        for iter = 1:Options.MaxIterations
          nu = mwFreq + gamman/h*B_;
          q = 1 - (aiso./(2*h*nu)).^2;
          Bnew = aiso./(gammae+gamman)./q.* ...
            (-mI+sign(aiso)*sqrt(mI.^2+q.*((h*nu/aiso).^2 - (I+1/2)^2)));
          if (B_~=0), RelChangeB = 1 - Bnew./B_; end
          B_ = Bnew;
          gamman = gn(iNucGrp)*nmagn; % re-include NZI
          if abs(RelChangeB)<Options.Accuracy, break; end
        end
        if (iter>maxIterationsDone), maxIterationsDone = iter; end
        if (iter==Options.MaxIterations)
          error(sprintf('Breit-Rabi solver didn''t converge after %d iterations!',iter));
        end
        Positions = [Positions B_];
        Intensities = [Intensities nn(iF)/nLines*ones(size(B_))];
      end
      logmsg(1,'  maximum %d iterations done',maxIterationsDone);
      
    else
      logmsg(1,'  perturbation expansion, order %d, accuracy %g',PerturbOrder,Options.Accuracy);
      
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
          for iter = 1:Options.MaxIterations
            if (abs(dB/Bb)<Options.Accuracy), break; end
            dB = polyval(c,Bb)/polyval(dc,Bb);
            Bb = Bb - dB;
          end
          if (iter>maxIterationsDone), maxIterationsDone = iter; end
          if (iter==Options.MaxIterations)
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
    if (PerturbOrder==0)
      logmsg(1,'  Breit-Rabi formula');
      
      % Breit-Rabi formula, J. A. Weil, J. Magn. Reson. 4, 394-399 (1971), [1]
      B_ = Exp.Field*1e-3; % T
      gebB = giso*bmagn*B_;
      gnbB = gn(iNucGrp)*nmagn*B_;
      for iF = 1:nFSpins
        if (nn(iF)==0), continue; end
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
          c(1) = giso*bmagn*B_;
          c(2) = aiso*mI;
          c(3) = pre(1) *      (I*(I+1)-mI^2);
          c(4) = pre(2) * mI * (1-2*I*(I+1)+2*mI^2);
          c(5) = pre(3) *      (I-2*I^3-I^4+6*(I^2+I-1)*mI^2-5*mI^4);
          c(6) = pre(4) * mI * (1+6*(I-1)*I*(I+1)*(I+2)-10*(-3+2*I*(I+1))*mI^2+14*mI^4);
          dE = sum(c(1:PerturbOrder+1));
          nu_ = dE/planck; % J -> Hz
          Positions(end+1) = nu_; % Hz
          Intensities(end+1) = nn(iF)/nLines;
        end
      end
      
    end

  end
  Shifts{p} = Positions - CentralResonance; % mT or Hz
  Amp{p} = Intensities;
  
  logmsg(1,'  spin group %d: %d F spins, %d lines',iNucGrp,nFSpins,numel(Shifts{p}));
  p = p + 1;
end

% Statistics: minimum, maximum, number of lines
%------------------------------------------------
nPeaks = 1;
posmax = CentralResonance;
posmin = CentralResonance;
for iShift = 1:length(Shifts)
  nPeaks = nPeaks * numel(Shifts{iShift});
  posmax = posmax + max(Shifts{iShift});
  posmin = posmin + min(Shifts{iShift});
end
logmsg(1,'  total %d lines',nPeaks);


if (FastMotionRegime)
  % Add Lorentzian
  LorentzianLw = FastMotionLw + Sys.lw(2);
  Sys.lw(2) = 0;
  maxLw = max(max(FastMotionLw),Sys.lw(1));
else
  maxLw = max(Sys.lw);
end

% Autoranging
%------------------------------------------------
if (SweepAutoRange)
  if FieldSweep
    posrange = (posmax-posmin)*Options.Stretch;
    Exp.Range = [posmin,posmax] + [-1 1]*max(5*maxLw,posrange);
    Exp.Range(1) = max(Exp.Range(1),0);
    logmsg(1,'  automatic field range from %g mT to %g mT',Exp.Range(1),Exp.Range(2));
  else
    posrange = (posmax-posmin)*Options.Stretch; % Hz
    Exp.Range = [posmin,posmax] + [-1 1]*max(5*maxLw/1e3,posrange);
    Exp.Range(1) = max(Exp.Range(1),0);
    Exp.Range = Exp.Range/1e9; % Hz -> GHz
    logmsg(1,'  automatic frequency range from %g GHz to %g GHz',Exp.Range(1),Exp.Range(2));
  end
end

if FieldSweep
  logmsg(1,'  spectral spread %g mT\n  g=%g resonance at %g mT',(posmax-posmin),giso,CentralResonance);
else
  logmsg(1,'  spectral spread %g MHz\n  g=%g resonance at %g GHz',(posmax-posmin)/1e6,giso,CentralResonance/1e9);
end

%===================================================================
% Spectrum construction
%===================================================================
logmsg(1,'Combining resonance shifts...');
if (nPeaks>1)
  Positions = sum(allcombinations(Shifts{:}),2) + CentralResonance;
  Intensity = prod(allcombinations(Amp{:}),2);
else
  Positions = CentralResonance;
  Intensity = 1;
end
if ~FieldSweep
  Positions = Positions/1e9; % Hz -> GHz
end

% Parallel mode: no intensities
if ParallelMode
  Intensity = Intensity*0;
end

% Temperature: include Boltzmann equilibrium polarization
if isfinite(Exp.Temperature)
  if FieldSweep
    e = exp(-planck*Exp.mwFreq*1e9/boltzm/Exp.Temperature);
    Population = [1 e]/(1+e);
    Polarization = Population(1)-Population(2);
    Intensity = Intensity*Polarization;
  else
    error('Boltzmann populations for frequency sweeps not implemented.');
  end
end

% Transition amplitudes
%--------------------------------------------------------------
g1mean2 = giso^2;
TransitionRate = (8*pi^2)*g1mean2*(bmagn/planck/1e9/2)^2;
Intensity = Intensity*TransitionRate;
if FieldSweep
  % 1/g factor (mT/MHz)
  dBdE = planck/(giso*bmagn)*1e9;
  Intensity = Intensity*dBdE;
end
%--------------------------------------------------------------

Harmonic2Do = Exp.Harmonic;

SweepRange = Exp.Range;
xAxis = linspace(SweepRange(1),SweepRange(2),Exp.nPoints);
Exp.deltaX = xAxis(2)-xAxis(1);

% (1) Fast-motion broadening
%---------------------------------------------------------
if (FastMotionRegime)
  logmsg(1,'Constructing spectrum with fast-motion Lorentzian linewidths...');

  dxx = min(Exp.deltaX,min(LorentzianLw)/5);
  nnPoints = round((SweepRange(2)-SweepRange(1))/dxx+1);
  xx = linspace(SweepRange(1),SweepRange(2),nnPoints);

  spec = 0;
  for iLine = 1:numel(Positions)
    spec = spec + Intensity(iLine)*lorentzian(xx,Positions(iLine),LorentzianLw(iLine),Harmonic2Do,Exp.mwPhase);
  end
  Exp.mwPhase = 0;
  Harmonic2Do = 0;
    
else
  
  xx = xAxis;
  dxx = Exp.deltaX;
  
  logmsg(1,'Constructing stick spectrum...'); 
  spec = constructspectrum(Positions,Intensity,SweepRange,Exp.nPoints);
  
end

spec = spec/Exp.deltaX;

% (2) Convolutional broadening
%---------------------------------------------------------
if (ConvolutionBroadening)
  
  fwhmL = Sys.lw(2);
  fwhmG = Sys.lw(1);
  if (fwhmL>0)
    HarmonicL = Harmonic2Do;
    mwPhaseL = Exp.mwPhase;
    HarmonicG = 0;
    mwPhaseG = 0;
  else
    HarmonicL = 0;
    mwPhaseL = 0;
    HarmonicG = Harmonic2Do;
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
  if (fwhmL>0)
    logmsg(1,'Convoluting with Lorentzian (FWHM %g %s, derivative %d)...',fwhmL,unitstr,HarmonicL);
    spec = convspec(spec,dxx,fwhmL,HarmonicL,0,mwPhaseL);
  end
  % Convolution with Gaussian
  if (fwhmG>0)
    logmsg(1,'Convoluting with Gaussian (FWHM %g %s, derivative %d)...',fwhmG,unitstr,HarmonicG);
    spec = convspec(spec,dxx,fwhmG,HarmonicG,1,mwPhaseG);
  end
  
end

if numel(xx)~=numel(xAxis)
  logmsg(1,'Re-interpolation (%d -> %d points)...',numel(xx),numel(xAxis));
  spec = interp1(xx,spec,xAxis);
end

% (3) Field modulation
%-----------------------------------------------------------------------
if FieldSweep
  if (Exp.ModAmp>0)
    logmsg(1,'  applying field modulation');
    spec = fieldmod(xAxis,spec,Exp.ModAmp);
  else
    % derivatives already included in the convolution
  end
else
  % frequency sweeps: modulation not available
end

%===================================================================

switch (nargout)
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
function Combs = allcombinations(varargin)

if (nargin==0), Combs = []; return; end

Combs = varargin{1}(:);
nCombs = numel(Combs);

for iArg = 2:nargin
  New = varargin{iArg}(:);
  nNew = numel(New);
  [idxNew,idxCombs] = find(ones(nNew,nCombs));
  Combs = [Combs(idxCombs(:),:), New(idxNew(:))];
  nCombs = nCombs*nNew;
end

return

%=========================================================
function Spectrum = constructspectrum(Positions,Amplitudes,Range,nPoints)

% Convert Positions to indices into spectral vector
idxPositions = (Positions-Range(1))/diff(Range) * (nPoints-1);
idxPositions = 1 + round(idxPositions);

% Remove out-of-range points
OutOfRange = (idxPositions<1) | (idxPositions>nPoints);
idxPositions(OutOfRange) = [];
Amplitudes(OutOfRange) = [];

% Bin all remaining lines into spectrum
Spectrum = full(sparse(1,idxPositions,Amplitudes,1,nPoints));

return
