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
FieldAutoRange = (~isfield(Exp,'Range') || isempty(Exp.Range)) && ...
  (~isfield(Exp,'CenterSweep') || isempty(Exp.CenterSweep));
if ~isfield(Options,'IsoCutoff'), Options.IsoCutoff = 1e-6; end

if ~isfield(Sys,'singleiso') || (Sys.singleiso==0)
  
  [SysList,weight] = expandcomponents(Sys,Options.IsoCutoff);
    
  if (numel(SysList)>1) && FieldAutoRange
    error('Multiple components: Please specify magnetic field range manually using Exp.Range or Exp.CenterSweep.');
  end
  
  spec = 0;
  for iComponent = 1:numel(SysList)
    [fieldAxis,spec_,Bres] = garlic(SysList{iComponent},Exp,Options);
    spec = spec + spec_*weight(iComponent);
  end
  
  % Output and plotting
  switch nargout
    case 0
      cla
      if (fieldAxis(2)<10000)
        plot(fieldAxis,spec);
        xlabel('magnetic field (mT)');
      else
        plot(fieldAxis/1e3,spec);
        xlabel('magnetic field (T)');
      end
      axis tight
      ylabel('intensity (arb.u.)');
      if isfield(Sys,'tcorr')
        fmStr = sprintf(', tcorr %g ns',Sys.tcorr*1e9);
      else
        fmStr = '';
      end
      title(sprintf('%0.8g GHz, %d points%s',...
        Exp.mwFreq,numel(fieldAxis),fmStr));
    case 1, varargout = {spec};
    case 2, varargout = {fieldAxis,spec};
    case 3, varargout = {fieldAxis,spec,Bres};
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
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% Checking experiment structure
%-------------------------------------------------------------------------
DefaultExp.Harmonic = 1;
DefaultExp.nPoints = 1024;
DefaultExp.ModAmp = 0;
DefaultExp.Mode = 'perpendicular';
DefaultExp.mwPhase = 0;
DefaultExp.Temperature = NaN; % don't compute thermal equilibrium polarizations

Exp = adddefaults(Exp,DefaultExp);

% Microwave frequency
if ~isfield(Exp,'mwFreq')
  error('Please supply the microwave frequency in Exp.mwFreq.');
end
if (numel(Exp.mwFreq)~=1) || any(Exp.mwFreq<=0) || ~isreal(Exp.mwFreq)
  error('Unintelligible microwave frequency in Exp.mwFreq.');
end

% Magnetic field range
if isfield(Exp,'CenterSweep')
  if isfield(Exp,'Range')
    logmsg(0,'Using Exp.CenterSweep and ignoring Exp.Range.');
  end
  Exp.Range = Exp.CenterSweep(1) + [-1 1]*Exp.CenterSweep(2)/2;
end

FieldAutoRange = ~isfield(Exp,'Range');
if ~FieldAutoRange
  if (Exp.Range(1)>=Exp.Range(2)) || any(Exp.Range<0)
    error('Invalid magnetic field range!');
  end
end

% Number of points
if any(~isreal(Exp.nPoints)) || numel(Exp.nPoints)>1 || (Exp.nPoints<2)
  error('Problem with Exp.nPoints. Needs to be a number not smaller than 2.')
end


% Detection harmonic
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
  logmsg(1,'  field modulation, amplitude %g mT',Exp.ModAmp);
  if (Exp.Harmonic<1)
    error('With field modulation (Exp.ModAmp), Exp.Harmonic=0 does not work.');
  end
  Exp.Harmonic = Exp.Harmonic - 1;
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

mwFreq = Exp.mwFreq*1e9; % GHz -> Hz
if Sys.fullg
  giso = mean(mean(Sys.g(1:3,:)));
else
  giso = mean(Sys.g(1,:));
end
CentralResonance = (mwFreq*planck)/(bmagn*giso)*1e3; % in the absence of nuclei,  mT

if (FastMotionRegime)
  FastMotionLw = fastmotion(Sys,CentralResonance,Sys.tcorr);
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
  % Make sure zero field splitting is smaller than mw frequency.
  if sum(abs(a_all).*I_all)*2>mwFreq
    disp('Microwave frequency is smaller than largest splitting of levels at zero field. Spectrum might be inaccurate. Consider using the function pepper instead.');
  end
  if any(I_all==0)
    error('Nuclei with spin 0 present.');
  end
  if any(a_all==0)
    error('Nuclei with coupling 0 present.');
  end
end

logmsg(1,'Computing resonance field shifts...');
p = 1;
Shifts = {};
Amp = {};

for iNucGrp = 1:Sys.nNuclei
  if (n_all(iNucGrp)==0), continue; end

  [I_this,nn] = equivcouple(I_all(iNucGrp),n_all(iNucGrp));
  nLines = (2*I_all(iNucGrp)+1)^n_all(iNucGrp);
  nFSpins = length(nn);
  a = a_all(iNucGrp)*planck;
  Positions = [];
  Intensities = [];

  if PerturbOrder==0
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
      B = 0;
      RelChangeB = inf;
      gamman = 0; % to obtain start field, neglect NZI
      for iter = 1:Options.MaxIterations
        nu = mwFreq + gamman/h*B;
        q = 1 - (a./(2*h*nu)).^2;
        Bnew = a./(gammae+gamman)./q.* ...
          (-mI+sign(a)*sqrt(mI.^2+q.*((h*nu/a).^2 - (I+1/2)^2)));
        if (B~=0), RelChangeB = 1 - Bnew./B; end
        B = Bnew;
        gamman = gn(iNucGrp)*nmagn; % re-include NZI
        if abs(RelChangeB)<Options.Accuracy, break; end
      end
      if (iter>maxIterationsDone), maxIterationsDone = iter; end
      if (iter==Options.MaxIterations)
        error(sprintf('Breit-Rabi solver didn''t converge after %d iterations!',iter));
      end
      Positions = [Positions B];
      Intensities = [Intensities nn(iF)/nLines*ones(size(B))];
    end
    logmsg(1,'  maximum %d iterations done',maxIterationsDone);
    
  else
    logmsg(1,'  perturbation expansion, order %d, accuracy %g',PerturbOrder,Options.Accuracy);
    
    % Based on \Delta E - h nu = 0 with a 5th-order Taylor expansion in a of
    % the Breit-Rabi expression for \Delta E. The resulting polynomial in a with terms
    % (0,a,a^2,a^3,a^4,...) is also a polynomial in B with the same terms
    % (B,0,1/B,1/B^2,1/B^3,...)
    maxIterationsDone = -1;
    pre = a .* (a/(giso*bmagn+gn(iNucGrp)*nmagn)/2).^(1:7);
    for iF = 1:nFSpins
      if nn(iF)==0, continue; end
      I = I_this(iF);
      for mI = -I:I
        c(6) = pre(4) * mI*(1+6*(I-1)*I*(I+1)*(I+2)-10*(-3+2*I*(I+1))*mI^2+14*mI^4);
        c(5) = pre(3) * (I-2*I^3-I^4 + 6*(I^2+I-1)*mI^2-5*mI^4);
        c(4) = pre(2) * mI*(1-2*I*(I+1)+2*mI^2);
        c(3) = pre(1) * (I*(I+1)-mI^2);
        c(2) = a*mI - mwFreq*planck;
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
        Positions(end+1) = Bb;
        Intensities(end+1) = nn(iF)/nLines;
      end
    end
    logmsg(1,'  maximum %d iterations done',maxIterationsDone);

  end
  
  Shifts{p} = Positions*1e3 - CentralResonance; % in mT
  Amp{p} = Intensities;
  p = p + 1;
  logmsg(1,'  spin group %d: %d F spins, %d lines',iNucGrp,nFSpins,numel(Positions));
end

% Statistics: minimum, maximum, number of lines
%------------------------------------------------
nPeaks = 1;
Bmax = CentralResonance;
Bmin = CentralResonance;
for k=1:length(Shifts)
  nPeaks = nPeaks * numel(Shifts{k});
  Bmax = Bmax + max(Shifts{k});
  Bmin = Bmin + min(Shifts{k});
end
logmsg(1,'  total %d lines',nPeaks);
%------------------------------------------------

if (FastMotionRegime)
  % Add Lorentzian
  LorentzianLw = FastMotionLw + Sys.lw(2);
  Sys.lw(2) = 0;
  maxLw = max(max(FastMotionLw),Sys.lw(1));
else
  maxLw = max(Sys.lw);
end

if (FieldAutoRange)
  Brange = (Bmax-Bmin)*Options.Stretch;
  Exp.Range = [Bmin,Bmax] + [-1 1]*max(5*maxLw,Brange);
  Exp.Range(1) = max(Exp.Range(1),0);
  logmsg(1,'  automatic field range from %g mT to %g mT',Exp.Range(1),Exp.Range(2));
end

logmsg(1,'  spectral spread %g mT\n  g=%g resonance at %g mT',(Bmax-Bmin),giso,CentralResonance);

%===================================================================
% Spectrum construction
%===================================================================
logmsg(1,'Combining resonance field shifts...');
if (nPeaks>1)
  B = sum(allcombinations(Shifts{:}),2) + CentralResonance;
  A = prod(allcombinations(Amp{:}),2);
else
  B = CentralResonance;
  A = 1;
end

% Parallel mode: no intensities
if ParallelMode
  A = A*0;
end

% Temperature: include Boltzmann equilibrium polarization
if isfinite(Exp.Temperature)
  e = exp(-planck*Exp.mwFreq*1e9/boltzm/Exp.Temperature);
  Population = [1 e]/(1+e);
  Polarization = Population(1)-Population(2);
  A = A*Polarization;
end

% Transition amplitudes
%--------------------------------------------------------------
% approximate only, but fairly accurate
% (correct expression would be the average over theta and phi of
% Eq. (A5) in Lund, Spectrochim Acta A 69 1294-1300 2008)
g1mean2 = giso^2;
TransitionRate = g1mean2*(bmagn/planck/1e9)^2;
% 1/g factor
dBdE = planck/(giso*bmagn)*1e9;
A = A*TransitionRate*dBdE*2*pi^2;
%--------------------------------------------------------------

Harmonic2Do = Exp.Harmonic;

FieldRange = Exp.Range;
x = linspace(FieldRange(1),FieldRange(2),Exp.nPoints);
Exp.deltaX = x(2)-x(1);

% (1) Fast-motion broadening
%---------------------------------------------------------
if (FastMotionRegime)
  logmsg(1,'Constructing spectrum with fast-motion Lorentzian linewidths...');

  dxx = min(Exp.deltaX,min(LorentzianLw)/5);
  nnPoints = round((FieldRange(2)-FieldRange(1))/dxx+1);
  xx = linspace(FieldRange(1),FieldRange(2),nnPoints);

  spec = 0;
  for iLine = 1:numel(B)
    spec = spec + A(iLine)*lorentzian(xx,B(iLine),LorentzianLw(iLine),Harmonic2Do);
  end
  Harmonic2Do = 0;
    
else
  
  xx = x;
  dxx = Exp.deltaX;
  
  logmsg(1,'Constructing stick spectrum...'); 
  spec = constructspectrum(B,A,FieldRange,Exp.nPoints);
  
end

spec = spec/Exp.deltaX;

% (2) Convolutional broadening
%---------------------------------------------------------
GaussianFWHM = Sys.lw(1);
if (GaussianFWHM>0)
  logmsg(1,'Convoluting with Gaussian (FWHM %g mT)...',GaussianFWHM);
  spec = convspec(spec,dxx,GaussianFWHM,Harmonic2Do,1);
  Harmonic2Do = 0;
end

LorentzianFWHM = Sys.lw(2);
if (LorentzianFWHM>0)
  logmsg(1,'Convoluting with Lorentzian (FWHM %g mT)...',LorentzianFWHM);
  spec = convspec(spec,dxx,LorentzianFWHM,Harmonic2Do,0,Exp.mwPhase);
  %Harmonic2Do = 0;
else
  if (Exp.mwPhase~=0)
    error('For Exp.mwPhase different from zero, please specify a Lorentzian broadening in Sys.lwpp or Sys.lw.');
  end
end

if numel(xx)~=numel(x)
  spec = interp1(xx,spec,x);
end

% (3) Field modulation
%-----------------------------------------------------------------------
if (Exp.ModAmp>0)
  logmsg(1,'  applying field modulation');
  spec = fieldmod(x,spec,Exp.ModAmp);
else
  % derivatives already included in the convolution
end

%===================================================================

switch (nargout)
  case 1, varargout = {spec};
  case 2, varargout = {x,spec};
  case 3, varargout = {x,spec,B};
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
