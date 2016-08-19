% saffron2    Simulate pulse EPR spectra
%
%     [x,S] = saffron(Sys,Exp,Opt)
%     [x,S,out] = saffron(Sys,Exp,Opt)
%
%     Sys   ... spin system with electron spin and ESEEM nuclei
%     Exp   ... experimental parameters (time unit us)
%     Opt   ... simulation options
%
%     out:
%       x       ... time or frequency axis (x{1} to x{n} for nD experiments)
%       S       ... simulated signal (ESEEM) or spectrum (ENDOR)
%       out     ... structure with FFT of ESEEM signal

function varargout = saffron2_full(Sys,Exp,Opt)

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
Link = 'epr@eth'; eschecker; error(LicErr); clear Link LicErr
% --------License ------------------------------------------------

% Guard against wrong number of input or output arguments.
if (nargin<2) || (nargin>3), error('Wrong number of input arguments!'); end
if (nargout<0), error('Not enough output arguments.'); end
if (nargout>4), error('Too many output arguments.'); end

% Initialize options structure to zero if not given.
if (nargin<3), Opt = struct('unused',NaN); end
if isempty(Opt), Opt = struct('unused',NaN); end

if ~isstruct(Exp)
  error('Second input argument (Exp) must be a structure!');
end
if ~isstruct(Opt)
  error('Third input argument (Opt) must be a structure!');
end

% User defined primary output for time dependent simulations
if ~isfield(Opt,'Output'), Opt.Output = 'Time'; end

if ~isfield(Opt,'Ordinate'), Opt.Ordinate = 'Complex'; end

% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

% Loop over species and isotopologues
%==================================================================

if ~isfield(Opt,'IsoCutoff'), Opt.IsoCutoff = 1e-4; end

if ~isfield(Sys,'singleiso') || (Sys.singleiso==0)
  
  % parse options
  [Opt.Output,err] = parseoption(Opt,'Output',{'Time','Frequency'});
  error(err);
  
  [Opt.Ordinate,err] = parseoption(Opt,'Ordinate',{'Complex','Real','Absolute'});
  error(err);
  
  [SysList,weight] = expandcomponents(Sys,Opt.IsoCutoff);
  
  ysum = 0;
  zsum = 0;
  for iComponent = 1:numel(SysList)
    [x,y_,out] = saffron2_full(SysList{iComponent},Exp,Opt);
    ysum = ysum + y_*weight(iComponent);
    zsum = zsum + out.fd*weight(iComponent);
  end
  
  out.td = ysum;
  switch Opt.Ordinate
    case 1, % complex output
    case 2, zsum = real(zsum);
    case 3, zsum = abs(zsum);
  end
  out.fd = zsum;
  
  switch Opt.Output
    case 1
      switch nargout
        case 0, % plotting, done below
        case 1, varargout = {ysum};
        case 2, varargout = {x,ysum};
        case 3, varargout = {x,ysum,out};
      end
    case 2
      switch nargout
        case 0, % plotting, done below
        case 1, varargout = {zsum};
        case 2, varargout = {out.f,zsum};
        case 3, varargout = {out.f,zsum,out};
      end
  end
  %===============================================================
  % Plotting
  %===============================================================
  ShowPlots = (nargout==0);
  if ShowPlots
    logmsg(1,'Graphical rendering...');
    clf
    
    if ~iscell(x)
      
      % Time domain
      subplot(2,1,1);
      plotQuadratureSignal = ~isreal(out.td);
      if plotQuadratureSignal
        h = plot(x,real(out.td),'b',x,imag(out.td),'r');
        set(h(1),'Color',[0 0 1]);
        set(h(2),'Color',[0.8 0.6 1]);
        legend('Re','Im');
        legend boxoff
      else
        plot(x,real(out.td));
      end
      axis tight
      xl = xlim;
      xlim([0 xl(2)]);
      yl = ylim;
      ylim(yl+[-1 1]*diff(yl)*0.1);
      
      xlabel('t (\mus)');
      ylabel('echo amplitude (arb.u.)');
      title('User-defined experiment, TD signal');
      set(gca,'Layer','top');
      
      % Frequency domain
      subplot(2,1,2);
      idx = find(out.f==0):length(out.f);
      xf = out.f(idx);
      if plotQuadratureSignal
        h = plot(xf,abs(out.fd(idx)),'g',xf,real(out.fd(idx)),'b',xf,imag(out.fd(idx)),'r');
        legend('abs','Re','Im');
        legend boxoff
        set(h(2),'Color',[0   0   1]);
        set(h(3),'Color',[0.8 0.6 1]);
        set(h(1),'Color',[0   0.8    0]);
      else
        h = plot(xf,abs(out.fd(idx)),'b',xf,real(out.fd(idx)),'g');
        set(h(2),'Color',[0 0.8 0]);
        legend('abs','Re');
        legend boxoff
      end
      axis tight
      xlim([0 max(out.f)]);
      xlabel('\nu (MHz)');
      ylabel('intensity (arb.u.)');
      title('Spectrum');
      
    elseif (iscell(x) && numel(x)==2)
      
      subplot(1,2,1);
      
      pcolor(x{1},x{2},real(out.td.')); shading flat; axis equal tight;
      set(gca,'Layer','top');
      title('Time domain (real part)');
      xlabel('t_1 (\mus)');
      ylabel('t_2 (\mus)');
      
      subplot(1,2,2);
      
      fx1 = fdaxis(Exp.dt(1),size(out.fd,1));
      if numel(Exp.dt)<2, Exp.dt(2) = Exp.dt(1); end
      fx2 = fdaxis(Exp.dt(2),size(out.fd,2));
      fd = abs(out.fd);
      if isfield(Opt,'logplot') && Opt.logplot
        fd = log(fd);
        maxfd = max(max(fd));
        fd(fd<maxfd-6) = maxfd-6;
      end
      pcolor(fx1,fx2,fd.'); shading flat; axis equal tight
      set(gca,'Layer','top');
      title('Frequency domain');
      xlabel('\nu_1 (MHz)');
      ylabel('\nu_2 (MHz)');
      fm1 = max(abs(fx1));
      fm2 = max(abs(fx2));
      fm = max(fm1,fm2);
      line([-1 1]*fm,[0 0],'Color','w');
      line([0 0],[-1 1]*fm,'Color','w');
      line([-1 1]*fm,[-1 1]*fm,'Color','w','LineStyle',':');
      line([1 -1]*fm,[-1 1]*fm,'Color','w','LineStyle',':');
      ylim([0 1]*fm);
    end
    
  end
  
  return
  
end
%==================================================================


logmsg(1,['=begin=saffron====' datestr(now) '=================']);
logmsg(2,'  log level %d',EasySpinLogLevel);
logmsg(1,'-general-----------------------------------------------');


%===================================================================
% Spin system
%===================================================================
if ~isfield(Sys,'Nucs'), Sys.Nucs = ''; end
outi = isotopologues(Sys.Nucs);
if outi.nIso>1
  error('saffron does not support isotope mixtures. Please specify pure isotopes in Sys.Nucs.');
end

[Sys,err] = validatespinsys(Sys);
error(err);
logmsg(1,'spins: %d electrons, %d nuclei',Sys.nElectrons,Sys.nNuclei);

maxNuclei = 40;
if (Sys.nNuclei>maxNuclei)
  error('saffron does not support systems with more than %d nuclei.',maxNuclei);
end

%===================================================================
% Experiment structure
%===================================================================

% Default settings
%===================================================================
DefaultExp.Temperature = [];
DefaultExp.Ordering = [];

DefaultExp.CrystalOrientation = [];
DefaultExp.CrystalSymmetry = '';
DefaultExp.MolFrame = [];

Exp = adddefaults(Exp,DefaultExp);

% Field
if ~isfield(Exp,'Field')
  error('Exp.Field is missing. Give a magnetic field in mT.');
end

% Temperature
if ~isempty(Exp.Temperature)
  error('Exp.Temperature is not supported for pulse EPR simulations.');
end

% Powder vs. crystal simulation
if isfield(Exp,'Orientation') || isfield(Exp,'Orientations')
  error('Exp.Orientation and Exp.Orientations are obsolete (as of EasySpin 5), use Exp.CrystalOrientation instead.');
end
PowderSimulation = isempty(Exp.CrystalOrientation);
Exp.PowderSimulation = PowderSimulation;

% Partial ordering
if ~isempty(Exp.Ordering)
  if ~PowderSimulation
    error('Partial ordering (Exp.Ordering) can only be used in a powder simulation.');
  else
    error('Partial ordering (Exp.Ordering) is not implemented in saffron.');
  end
end

% T1, T2
if ~isfield(Sys,'T1T2')
  if isfield(Exp,'T1T2')
    error('T1 and T2 should be provided in Sys.T1T2 and not in Exp.T1T2.');
  end
  Sys.T1T2 = [0 0];
end
Sys.T1T2(Sys.T1T2==0) = inf;
if numel(Sys.T1T2)~=2
  error('Sys.T1T2 must contain two numbers, T1 and T2 in microseconds.');
end
if any(Sys.T1T2<=0) || any(~isreal(Sys.T1T2))
  error('T1 and T2 in Sys.T1T2 must be positive, in microseconds.');
end

% Pulse sequence
%===================================================================
% User-specified pulse sequence -----------------------------------------
logmsg(1,'User-specified pulse experiment.');

if isfield(Exp,'Sequence')
  error('Predefined pulse sequences not supported in saffron2_full.');
end

if any(~isinf(Sys.T1T2))
  error('Relaxation times T1 and T2 for custom sequences not supported.');
end

% Pulse parameters
if ~isfield(Exp,'Flip') && ~isfield(Exp,'Amplitude') && ~isfield(Exp,'Pulse')
  error('For a pulse experiment, give either Exp.Sequence or define the pulse sequence explicitly.');
end
if isfield(Exp,'Flip')
  nIntervals = numel(Exp.Flip);
  if (~any(mod(Exp.Flip,1)) && max(Exp.Flip)<=4 && any(Exp.Flip))
    warning('The input to Exp.Flip is now in radians. Multiply the old Exp.Flip input by pi/2.')
  end
elseif isfield(Exp,'Pulse')
  nIntervals = numel(Exp.Pulse);
elseif isfield(Exp,'Amplitude')
  nIntervals = numel(Exp.tp);
else
  error('Exp.Flip or the Exp.Pulse structure needs to be defined for a custom pulse sequence.');
end

% Incrementation scheme
if ~isfield(Exp,'Inc')
  Exp.Inc(1:nIntervals) = 0; % echo detection if no incrementation scheme is given
end
if numel(Exp.Inc)~=nIntervals
  error('Exp.Inc must contain the same number of elements as Exp.Flip/Exp.Amplitude.');
end
IncScheme = Exp.Inc(Exp.Inc~=0); % ? include IncFreq/Inctp (field sweep separate)
if isempty(IncScheme), IncSchemeID = 0; % FID/echo detection for any type of pulse sequence
elseif isequal(IncScheme,1), IncSchemeID = 1;
elseif isequal(IncScheme,[1 1]), IncSchemeID = 2;
elseif isequal(IncScheme,[1 -1]), IncSchemeID = 3;
elseif isequal(IncScheme,[1 1 -1 -1]), IncSchemeID = 4;
elseif isequal(IncScheme,[1 -1 -1 1]), IncSchemeID = 5;
elseif isequal(IncScheme,[1 -1 1 -1]), IncSchemeID = 6;
elseif isequal(IncScheme,[1 2]), IncSchemeID = 11;
elseif isequal(IncScheme,[1 2 1]), IncSchemeID = 12;
elseif isequal(IncScheme,[1 2 2]), IncSchemeID = 13;
elseif isequal(IncScheme,[1 1 2]), IncSchemeID = 14;
elseif isequal(IncScheme,[1 2 2 1]), IncSchemeID = 15;
elseif isequal(IncScheme,[1 2 -2 1]), IncSchemeID = 16;
elseif isequal(IncScheme,[1 1 2 2]), IncSchemeID = 17;
elseif isequal(IncScheme,[1 -1 2]), IncSchemeID = 18;
elseif isequal(IncScheme,[1 1 -1 -1 2]), IncSchemeID = 19;
elseif isequal(IncScheme,[1 -1 -1 1 2]), IncSchemeID = 20;
elseif isequal(IncScheme,[1 -1 1 -1 2]), IncSchemeID = 21;
else
  error('Unsupported incrementation scheme in Exp.Inc.');
end
nDimensions = max(abs(Exp.Inc)); % ? check other Inc parameters

% Pulse delays
if ~isfield(Exp,'t')
  Exp.t = zeros(1,nIntervals);
end
if numel(Exp.t)~=nIntervals
  error('Exp.t must contain the same number of elements as Exp.Flip');
end
if all(Exp.t==0) && (numel(IncScheme)<nIntervals) % ? change, see above
  logmsg(0,'Some delays are zero, but are not incremented!');
end

% Detection
% ???????? detection operator
if ~isfield(Exp,'DetectionIntegrate')
  Exp.DetectionIntegrate = 0;
end
% Integration/detection window
if ~isfield(Exp,'DetectionStep')
  Exp.DetectionStep = 0.001; % us
end
if ~isfield(Exp,'DetectionDelay')
  Exp.DetectionDelay = 0; % us
end
if isfield(Exp,'DetectionWindow')
  Exp.DetectionPoints = round(Exp.DetectionWindow/Exp.DetectionStep);
elseif isfield(Exp,'DetectionPoints')
  Exp.DetectionWindow = Exp.DetectionPoints*Exp.DetectionStep;
end
if (Exp.DetectionIntegrate==1 && ~isfield(Exp,'DetectionWindow'))
  error('Echo integration is requested, but the integration window is not defined.')
end
if (isfield(Exp,'DetectionWindow') && (Exp.t(end)+Exp.DetectionDelay)<Exp.DetectionWindow/2)
  error('Detection window overlaps with the last pulse. Adjust window length or select the full transient option.');
end
% Options for detection
if (Exp.DetectionIntegrate==0 && ~isfield(Exp,'DetectionWindow'))
  Exp.Detection = 'singlepoint';
  logmsg(1,'Single point detection.');
elseif (Exp.DetectionIntegrate==1 && isfield(Exp,'DetectionWindow'))
  Exp.Detection = 'echodetection';
  logmsg(1,'Echo integration over a %d ns window in %d ns steps.',Exp.DetectionWindow*10^3,Exp.DetectionStep*10^3);
elseif (Exp.DetectionIntegrate==0 && isfield(Exp,'DetectionWindow'))
  Exp.Detection = 'echodetection';
  logmsg(1,'Echo transient acquisition in a %d ns window with %d ns time steps.',Exp.DetectionWindow*10^3,Exp.DetectionStep*10^3);
end
if ~isfield(Exp,'DetectionFullTransient')
  Exp.DetectionFullTransient = 0;
end
if Exp.DetectionFullTransient==1 % ? not implemented yet
  Exp.Detection = 'fulltransient';
  logmsg(1,'Transient detection of the full pulse sequence.');
end

if (isempty(IncScheme) && ...
    ~(strcmp(Exp.Detection,'echodetection') || strcmp(Exp.Detection,'fulltransient')))
  error('No incrementation scheme is given. For transient simulations specify an integration window.')
end

% dt
if ~isfield(Exp,'dt')
  if isempty(IncScheme) % echo detection
    Exp.dt = [];
  else
    error('Exp.dt is missing.');
  end
end
if numel(Exp.dt)==1
  Exp.dt = Exp.dt*ones(1,nDimensions);
elseif numel(Exp.dt)~=nDimensions
  error('Exp.dt needs either 1 or %d elements, one per dimension. You gave %d.',nDimensions,numel(Exp.dt));
end

% nPoints % ? check for pulse crossings
if ~isfield(Exp,'nPoints')
  if isENDOR
    Exp.nPoints = 1001;
  else
    if (nDimensions==0)
      Exp.nPoints = [];
    elseif (nDimensions==1)
      Exp.nPoints = 512;
    else
      Exp.nPoints = [1 1]*256;
    end
  end
end
if numel(Exp.nPoints)==1
  Exp.nPoints = Exp.nPoints*ones(1,nDimensions);
elseif numel(Exp.nPoints)~=nDimensions
  error('Exp.nPoints needs either 1 or %d elements, one per dimension. You gave %d.',nDimensions,numel(Exp.nPoints));
end

% Add detection delay to the last delay
if isfield(Exp,'t')
  Exp.t(end) = Exp.t(end)+Exp.DetectionDelay;
end
% Update settings if echo detection is used
if strcmp(Exp.Detection,'echodetection')
  nPoints = [Exp.nPoints Exp.DetectionPoints];
  dt = [Exp.dt Exp.DetectionStep];
  nDimensions = nDimensions+1; % Additional dimension for echo detection
  % Update incrementation scheme by adding additional final dimension
  if isempty(IncScheme)
    IncScheme = 1;
  else
    IncScheme = [IncScheme max(IncScheme)+1];
  end
  
  if isempty(IncScheme), IncSchemeID = 0; % FID/echo detection for any type of pulse sequence
  elseif isequal(IncScheme,1), IncSchemeID = 1;
  elseif isequal(IncScheme,[1 1]), IncSchemeID = 2;
  elseif isequal(IncScheme,[1 -1]), IncSchemeID = 3;
  elseif isequal(IncScheme,[1 1 -1 -1]), IncSchemeID = 4;
  elseif isequal(IncScheme,[1 -1 -1 1]), IncSchemeID = 5;
  elseif isequal(IncScheme,[1 -1 1 -1]), IncSchemeID = 6;
  elseif isequal(IncScheme,[1 2]), IncSchemeID = 11;
  elseif isequal(IncScheme,[1 2 1]), IncSchemeID = 12;
  elseif isequal(IncScheme,[1 2 2]), IncSchemeID = 13;
  elseif isequal(IncScheme,[1 1 2]), IncSchemeID = 14;
  elseif isequal(IncScheme,[1 2 2 1]), IncSchemeID = 15;
  elseif isequal(IncScheme,[1 2 -2 1]), IncSchemeID = 16;
  elseif isequal(IncScheme,[1 1 2 2]), IncSchemeID = 17;
  elseif isequal(IncScheme,[1 -1 2]), IncSchemeID = 18;
  elseif isequal(IncScheme,[1 1 -1 -1 2]), IncSchemeID = 19;
  elseif isequal(IncScheme,[1 -1 -1 1 2]), IncSchemeID = 20;
  elseif isequal(IncScheme,[1 -1 1 -1 2]), IncSchemeID = 21;
  else
    error('Echo detection/integration not supported for the current incrementation scheme.');
  end
  
else
  nPoints = Exp.nPoints;
  dt = Exp.dt;
end

% Pulse definitions
if isfield(Exp,'Phase') && (~any(mod(Exp.Phase,1)) && max(Exp.Phase)<=4 && any(Exp.Phase))
  warning('The input to Exp.Phase is now in radians. Multiply the old Exp.Phase input by pi/2.')
end
% Pulse lengths
if isfield(Exp,'tp') % Exp.tp input given, convert to Exp.Pulse
  
  if numel(Exp.tp)~=nIntervals
    error('Exp.tp contains a wrong number of elements.');
  end
  if ~isfield(Exp,'Pulse')
    Exp.Pulse = repmat(struct,nIntervals,1);
  end
  if numel(Exp.Pulse)<numel(Exp.tp)
    for i = numel(Exp.Pulse):numel(Exp.tp)
      Exp.Pulse(i).tp = Exp.tp(i);
    end
  end
  for p = 1:nIntervals
    if ~isfield(Exp.Pulse(p),'tp') || isempty(Exp.Pulse(p).tp)
      Exp.Pulse(p).tp = Exp.tp(p);
    elseif Exp.Pulse(p).tp~=Exp.tp(p)
      error('Pulse lengths in Exp.tp and Exp.Pulse(%i).tp do not agree.',p)
    end
  end
  if ~isfield(Exp,'Phase')
    Exp.Phase = (pi/2)*ones(1,nIntervals);  % y phase by default
  end
  
elseif ~isfield(Exp,'tp') && isfield(Exp,'Pulse') % pulse sequence only defined in Exp.Pulse
  
  if numel(Exp.Pulse)~=nIntervals
    error('Exp.Pulse contains a wrong number of elements.');
  end
  for p = 1:nIntervals
    Exp.tp(p) = Exp.Pulse(p).tp;
  end
  
else % ideal pulses
  
  Exp.tp = zeros(1,nIntervals);
  if ~isfield(Exp,'Phase')
    Exp.Phase = (pi/2)*ones(1,nIntervals);  % y phase by default
  end
end

% Determine whether to use ideal pulse theory
idealPulse = (Exp.tp==0);
realPulse = ~idealPulse;
if any(realPulse)
  for p = find((1:numel(Exp.tp)).*realPulse)
    
    if isfield(Exp,'Flip')
      if ~isfield(Exp.Pulse(p),'Flip') || isempty(Exp.Pulse(p).Flip)
        Exp.Pulse(p).Flip = Exp.Flip(p);
      elseif Exp.Pulse(p).Flip~=Exp.Flip(p)
        error('Flip angles in Exp.Flip and Exp.Pulse(%i).Flip do not agree.',p)
      end
    end
    
    if isfield(Exp,'Amplitude')
      if numel(Exp.Amplitude)==1
        Exp.Pulse(p).Amplitude = Exp.Amplitude;
      else
        Exp.Pulse(p).Amplitude = Exp.Amplitude(p);
      end
    end
    
    if isfield(Exp,'Phase')
      if ~isfield(Exp.Pulse(p),'Phase') || isempty(Exp.Pulse(p).Phase)
        Exp.Pulse(p).Phase = Exp.Phase(p);
      elseif Exp.Pulse(p).Phase~=Exp.Phase(p)
        error('Phase in Exp.Phase and Exp.Pulse(%i).Phase do not agree.',p)
      end
    else
      if ~isfield(Exp.Pulse(p),'Phase') || isempty(Exp.Pulse(p).Phase)
        Exp.Pulse(p).Phase = pi/2; % y phase by default
      end
    end
    
    % Time step for pulse propagation
    if isfield(Exp,'Pulsedt') && ~isempty(Exp.Pulsedt)
      if numel(Exp.Pulsedt)==nIntervals
        Exp.Pulse(p).TimeStep = Exp.Pulsedt(p);
      else
        Exp.Pulse(p).TimeStep = Exp.Pulsedt;
      end
    end
    
    [Exp.tpulse{p},Exp.IQpulse{p}] = pulse(Exp.Pulse(p));
    
  end
end

% Validate incrementation scheme: assert all delays are nonnegative
for iInterval = 1:nIntervals
  iDim = abs(Exp.Inc(iInterval));
  if (iDim==0)  % delay is kept constant
    t_range = Exp.t(iInterval)*[1 1];
  else % delay is incremented/decremented
    t_range = Exp.t(iInterval) + sign(Exp.Inc(iInterval))*(Exp.nPoints(iDim)-1)*Exp.dt(iDim);
  end
  if any(t_range<0)
    error('Negative delay after pulse %d.',iInterval);
  end
end

% Orientation selection
% Options: - orientation selection through simulation with real pulses
%          - ideal pulses, no orientation selection
if any(realPulse)
  if ~isfield(Exp,'mwFreq')
    error('Exp.mwFreq is required for simulations with real pulses. Please give Exp.mwFreq.');
  end
%   if isfield(Exp,'ExciteWidth')
%     warning('Exp.ExciteWidth is obsolete for simulations with real pulses. See documentation.')
%   end
end
if isfield(Exp,'ExciteWidth')
  warning('Exp.ExciteWidth is not used in saffron2_full.')
end
if any(realPulse)
  logmsg(1,'Orientation selection with real pulses.');
else
  logmsg(1,'No orientation selection (infinite bandwidth).');
end

if isfield(Exp,'HStrain')
  error('You gave Exp.HStrain, but it should be Sys.HStrain (in the system, not the experiment structure).');
end

%===================================================================
% Options structure
%===================================================================
%
if ~isfield(Opt,'Symmetry'), Opt.Symmetry = []; end
if ~isfield(Opt,'SymmFrame'), Opt.SymmFrame = []; end
if ~isfield(Opt,'Transitions'), Opt.Transitions = []; end
if ~isfield(Opt,'Sites'), Opt.Sites = []; end

% Expansion factor: determines size of spectral buffer
if ~isfield(Opt,'Expand')
  if (nDimensions==1), Opt.Expand = 4; else Opt.Expand = 2; end
end
maxExpand = 8;
if (numel(Opt.Expand)~=1) || (Opt.Expand<0) || (Opt.Expand>maxExpand) || rem(Opt.Expand,1)
  error('Opt.Expand must be an integer between 0 and %d.',maxExpand);
end
if ~isfield(Opt,'DetectionExpand')
  Opt.DetectionExpand = 4;
end

% Number of knots: determines number of orientations
if ~isfield(Opt,'nKnots'), Opt.nKnots = 30+1; end
if numel(Opt.nKnots)>1
  error('Only one number allowed in Opt.nKnots. saffron does not support interpolation.');
end
if (Opt.nKnots<7)
  error('Opt.nKnots must be at least 7. You gave %d.',Opt.nKnots);
end

if ~isfield(Opt,'OriThreshold'), Opt.OriThreshold = 0.005; end

if ~isfield(Opt,'Window'),
  if (nDimensions==1), Opt.Window = 'ham+'; else Opt.Window = 'ham+'; end
end

if ~isfield(Opt,'ZeroFillFactor');
  Opt.ZeroFillFactor = 2;
end

DataProcessing = 1;

% undocumented options
if ~isfield(Opt,'logplot'), Opt.logplot = 0; end
if ~isfield(Opt,'PartialIFFT'), Opt.PartialIFFT = 1; end
if ~isfield(Opt,'TimeDomain'), Opt.TimeDomain = 1; end


%==========================================================================
% Symmetry determination and orientational grid.
%==========================================================================
[Exp,Opt] = p_symandgrid(Sys,Exp,Opt);

% Process crystal orientations, crystal symmetry, and frame transforms
% This sets Orientations, nOrientations, nSites and AverageOverChi
[Orientations,nOrientations,nSites,AverageOverChi] = p_crystalorientations(Exp,Opt);

logmsg(1,'-Hamiltonians------------------------------------------');

%=====================================================================
% Compute spin Hamiltonian
%=====================================================================

logmsg(1,'setting up spin Hamiltonian...');
% coreSys = nucspinrmv(Sys,shfNuclei);
coreSys = Sys; % ????????
coreSys.processed = 0;
coreSys.lw = 0;
% Operators for constructing Hamiltonian
[F,Gx,Gy,Gz] = sham(coreSys);
% Operators for computing <i|S|i>
logmsg(1,'Generating all cartesian spin operators');
for iSpin = 1:numel(coreSys.Spins)
  SpinOps{iSpin,1} = sop(coreSys.Spins,iSpin,1);
  SpinOps{iSpin,2} = sop(coreSys.Spins,iSpin,2);
  SpinOps{iSpin,3} = sop(coreSys.Spins,iSpin,3);
end
% Total electron spin operators
for i = 1:3
  totSpinOps{i} = zeros(size(F));
end
for iSpin = 1:numel(coreSys.S)
  for i = 1:3
    totSpinOps{i} = totSpinOps{i} + SpinOps{iSpin,i};
  end
end

%=====================================================================

%=====================================================================
% Preparation
%=====================================================================
logmsg(1,'Preparation...');

if Opt.TimeDomain
  logmsg(1,'  time domain simulation');
else
  logmsg(1,'  frequency domain simulation');
  % Prepare for binning method
  if ~strcmp(Exp.Detection,'echodetection')
    ExpansionFactor = 2.^Opt.Expand;
    nPointsF = ExpansionFactor*nPoints;
    if numel(nPointsF)==1
      logmsg(1,'  %d points, expand x%d -> %d points',nPoints,ExpansionFactor,nPointsF);
    elseif numel(nPointsF)==2
      logmsg(1,'  %dx%d points, expand x%d -> %dx%d points',...
        nPoints(1),nPoints(2),ExpansionFactor,nPointsF(1),nPointsF(2));
    elseif numel(nPointsF)==3
      logmsg(1,'  %dx%dx%d points, expand x%d -> %dx%dx%d points',...
        nPoints(1),nPoints(2),nPoints(3),ExpansionFactor,nPointsF(1),nPointsF(2),nPointsF(3));
    end
  else % echo detection (increased expansion factor required for accurate echo simulation, at least 8)
    ExpansionFactor = 2.^Opt.Expand;
    nPointsF = ExpansionFactor*nPoints(1:numel(nPoints)-1);
    ExpansionFactorDet = (2.^Opt.DetectionExpand);
    nPointsF(numel(nPoints)) = ExpansionFactorDet*nPoints(end);
    if numel(nPointsF)==1
      logmsg(1,'  %d points, expand x%d -> %d points',nPoints,ExpansionFactorDet,nPointsF);
    elseif numel(nPointsF)==2
      logmsg(1,'  %dx%d points, expand x%d/x%d -> %dx%d points',...
        nPoints(1),nPoints(2),ExpansionFactor,ExpansionFactorDet,nPointsF(1),nPointsF(2));
    elseif numel(nPointsF)==3
      logmsg(1,'  %dx%dx%d points, expand x%d/x%d/x%d -> %dx%dx%d points',...
        nPoints(1),nPoints(2),nPoints(3),ExpansionFactor,ExpansionFactor,ExpansionFactorDet,nPointsF(1),nPointsF(2),nPointsF(3));
    end
  end
  
  % allocate array(s)
  if (nDimensions==1), siz = [1, nPointsF]; else siz = nPointsF; end
  buff = zeros(siz);
  buff(1) = 1e-300i; % make sure it's complex
  
end

if nDimensions==1
  totaltd = zeros(nPoints,1);
else
  totaltd = zeros(nPoints);
end
%=====================================================================


%=====================================================================
% Orientation loop
%=====================================================================
logmsg(1,'Looping over %d orientations...',nOrientations);

% Line broadening and offset integration
if any(realPulse)
  if ~isfield(Opt,'nOffsets');
    Opt.nOffsets = 21; % ? optimization loop
  end
  if isfield(Opt,'lwOffset')
    error('Opt.lwOffset is obsolete, the offset linewidth is derived from the HStrain parameter.');
  end
  if Opt.nOffsets==0 || Opt.nOffsets==1
    Opt.nOffsets = 1;
    offsets = zeros(nOrientations,1);
    offsetWeight= ones(nOrientations,1);
  else
    % Linewidths defined through HStrain for all orientations
    if isfield(Sys,'HStrain') && any(Sys.HStrain~=0)
      Sys.HStrain(Sys.HStrain==0) = 0.1;
      zLab = ang2vec(Orientations(:,1),Orientations(:,2));
      lwOffset = sum(diag(Sys.HStrain)*zLab,1);
      offsets = bsxfun(@times,repmat(linspace(-1,1,Opt.nOffsets),nOrientations,1),2*lwOffset.');
      offsetWeight = exp(-4*log(2)*(bsxfun(@rdivide,offsets,lwOffset.')).^2);
      offsetWeight = bsxfun(@rdivide,offsetWeight,sum(offsetWeight,2));
    else
      lw = 0.1; % MHz, default (????)
      offsets = repmat(linspace(-1,1,Opt.nOffsets)*lw*2,nOrientations,1);
      offsetWeight = exp(-4*log(2)*(offsets/lw).^2);
      offsetWeight = bsxfun(@rdivide,offsetWeight,sum(offsetWeight,2));
    end
  end
else
  Opt.nOffsets = 1;
  offsets = zeros(nOrientations,1);
  offsetWeight= ones(nOrientations,1);
end

nSkippedOrientations = 0;
for iOri = 1:nOrientations
  
  % Get magnetic field orientation
  %------------------------------------------------------------------
  [xLab_M,yLab_M,zLab_M] = erot(Orientations(iOri,:),'rows');
  % xLab_M, yLab_M, zLab_M represented in the molecular frame

  % total spin operators in the lab frame
  totSpinOps_L{1} = xLab_M(1)*totSpinOps{1} + xLab_M(2)*totSpinOps{2} + xLab_M(3)*totSpinOps{3};
  totSpinOps_L{2} = yLab_M(1)*totSpinOps{1} + yLab_M(2)*totSpinOps{2} + yLab_M(3)*totSpinOps{3};
  totSpinOps_L{3} = zLab_M(1)*totSpinOps{1} + zLab_M(2)*totSpinOps{2} + zLab_M(3)*totSpinOps{3};

  % Compute electronic Hamiltonian, energies and <S>
  %------------------------------------------------------------------
  H = F + Exp.Field*(zLab_M(1)*Gx + zLab_M(2)*Gy + zLab_M(3)*Gz);   
  [eV,eE] = eig((H+H')/2);
  
  % total spin operators in Hamiltonian eigenframe
  for i = 1:3
    totSpinOps_eig{i} = eV'*totSpinOps_L{i}*eV;
  end
  
  % ???? Detection operator
  Det = totSpinOps_eig{1} + 1i*totSpinOps_eig{2};
  
  % Equilibrium density matrix
  if ~isempty(Exp.Temperature)
    % make temperature aware, include Boltzmann populations
    Sig0 = sigeq(eE,Exp.Temperature);
  else
    % high temperature approximation?
    Sig0 = -totSpinOps_eig{3}; % -Sz
  end
  
  % Rotating frame transformation
  H0_rot = eE - Exp.mwFreq*1e3*totSpinOps_eig{3};
   
  logmsg(2,'orientation %d of %d',iOri,nOrientations);
   
  for iOffset = 1:Opt.nOffsets
    
    H_rot = H0_rot + offsets(iOri,iOffset)*totSpinOps_eig{3};
    E = diag(H_rot); % ! wrong???
    
    % Pulse propagators
    for iInt = 1:nIntervals
      if realPulse(iInt)
        % FullPulsePropagator = sf_propagator(taxis,IQpulse,H0,Sx,Sy)
        FullPulsePropagator{iInt} = sf_propagator(Exp.tpulse{iInt},Exp.IQpulse{iInt},H_rot,totSpinOps_eig{1},totSpinOps_eig{2});
      else % ideal pulses
        % FullPulsePropagator = expm(-1i*flipangle*(cos(phase)*Sx+sin(phase)*Sy))
        Op = cos(Exp.Phase(iInt))*totSpinOps_eig{1} + ...
             sin(Exp.Phase(iInt))*totSpinOps_eig{2};
        FullPulsePropagator{iInt} = expm(-1i*Exp.Flip(iInt)*Op);
      end
    end
    
    prefactor = Exp.OriWeights(iOri)*offsetWeight(iOri,iOffset);
        
    % Compute peaks and amplitudes
    % General method ------------------------------------------
    increments = Exp.Inc;
    eyeN = eye(size(Sig0));
    Block = cell(0);
    iBlock = 0;
    Left = prefactor*Sig0;
    Right = eyeN;
    for iInt = 1:nIntervals
      
      % Pulse propagator
      Left = FullPulsePropagator{iInt}*Left;
      Right = Right*FullPulsePropagator{iInt}';
      
      % Free evolution propagator
      if iInt==nIntervals && strcmp(Exp.Detection,'echodetection')
        t = Exp.t(iInt) - Exp.DetectionWindow/2;
      else
        t = Exp.t(iInt);
      end
      
      if (t>0)
        Left = diag(exp(-2i*pi*E*t)) * Left;
        Right = Right * diag(exp(+2i*pi*E*t));
%         Left = expm(-2i*pi*H_rot*t) * Left;
%         Right = Right * expm(+2i*pi*H_rot*t);
      end
      
      % Conclude block and start next one
      if increments(iInt)~=0
        if (iBlock==0)
          G = Left*Right;
          if numel(G)==1
            G = G*eyeN;
          end
        else
          if numel(Left)==1
            Block{iBlock} = Left*eyeN;
          else
            Block{iBlock} = Left;
          end
        end
        iBlock = iBlock + 1;
        Left = 1;
        Right = 1;
      end
      
    end % iInt
    
    % Additional interval for echo detection
    if strcmp(Exp.Detection,'echodetection')
      if (iBlock==0)
        G = Left*Right;
        if numel(G)==1
          G = G*eyeN;
        end
      else
        if numel(Left)==1
          Block{iBlock} = Left*eyeN;
        else
          Block{iBlock} = Left;
        end
      end
      Left = eyeN;
      Right = eyeN;
    end
    
    if increments(end)~=0
      D = Det;
    else
      D = Right*Det*Left;
    end
    
    % Accumulate peaks / generate time domain
    if Opt.TimeDomain
      
      totaltd = totaltd + ...
        sf_evolve_full(G,D,diag(E),nPoints,dt,IncScheme,Block);
%       totaltd = totaltd + ...
%         evolve(G,D,H_rot,nPoints,dt,IncScheme,Block);
      
    else
      
      if isempty(Block)
        BlockL = {};
        BlockR = {};
      else
        for iBlock = 1:numel(Block)
          BlockL{iBlock} = Block{iBlock};
          BlockR{iBlock} = Block{iBlock}';
        end
      end
      sf_peaks_full(IncSchemeID,buff,dt,E,G,D,BlockL{:},BlockR{:});
      
    end
    
  end % offset loop
  
end % orientation loop

logmsg(1,'end of orientation/transition loop');
logmsg(1,'%d of %d orientations skipped',nSkippedOrientations,nOrientations);
%=================================================================



%=================================================================
% Postprocessing
%=================================================================

if Opt.TimeDomain
  td = totaltd;
else
  logmsg(1,'Postprocessing...');
  if (nDimensions==1)
    td = ifft(buff)*numel(buff);
    td = td(1:nPoints).';
  elseif (nDimensions==2)
    td = ifft2dpartial(buff,nPoints,Opt.PartialIFFT)*numel(buff);
  end
end

% Normalize modulation signal
% No need to normalize out the experiment prefactors (due to
% pulse transfer amplitudes), since they were not included above.
EqDensityTrace = prod(2*Sys.I+1); %??????????
td = td/EqDensityTrace;
% td = td/nPathways;

% td = td/max(abs(pathwayprefactor));

% Echo integration
if Exp.DetectionIntegrate==1
  if nDimensions==1
    td = sum(td);
  else
    td = sum(td,nDimensions);
  end
  dt = dt(1:end-1);
  nPoints = nPoints(1:end-1);
  nDimensions = nDimensions-1;
end

%=================================================================
% Time axes, relaxation
%=================================================================
switch (nDimensions)
  case 0 % e.g. echo-detected field sweep (B in external loop)
    t = [];
  case 1
    t = (0:nPoints-1)*dt;
  case 2,
    clear t
    t{1} = (0:nPoints(1)-1)*dt(1);
    t{2} = (0:nPoints(2)-1)*dt(2);
  case 3,
    clear t
    t{1} = (0:nPoints(1)-1)*dt(1);
    t{2} = (0:nPoints(2)-1)*dt(2);
    t{3} = (0:nPoints(3)-1)*dt(3);
end

%===============================================================
% TD data processing
%===============================================================
logmsg(1,'-final-------------------------------------------------');
logmsg(1,'Data processing...');
switch nDimensions
  case 0
    fd = [];
    f = [];
  case 1
    if DataProcessing
      
      tdx = td;
      
      % Baseline correction
      tdx = tdx - mean(tdx);
      
      % Apodization
      win = apowin(Opt.Window,numel(tdx));
      tdx = tdx.*win;
      
      % Fourier transformation
      fd = fft(tdx,Opt.ZeroFillFactor*numel(tdx));
      fd = fftshift(fd);
      f = fdaxis(dt,length(fd));
    end
    
    
  case 2
    if DataProcessing
      
      tdx = basecorr(td,[1 2],[0 0]);
      
      w1 = apowin(Opt.Window,nPoints(1));
      w2 = apowin(Opt.Window,nPoints(2));
      
      fd = fftshift(fftn(tdx.*(w1*w2.'),Opt.ZeroFillFactor*nPoints));
      f{1} = fdaxis(dt(1),size(fd,1));
      f{2} = fdaxis(dt(2),size(fd,2));
    end
    
end

if max(abs(fd))<1e-300;
  fd = fd*0;
end

% Collect output structure
if DataProcessing
  out.f = f;
  out.fd = fd;
else
  out = [];
end


%===============================================================
% Output
%===============================================================
switch nargout
  case 1, varargout = {td};
  case 2, varargout = {t,td};
  case 3, varargout = {t,td,out};
end


%===============================================================
% Report performance
%===============================================================
[Hours,Minutes,Seconds] = elapsedtime(StartTime,clock);
if (Hours>0)
  msg = sprintf('saffron took %dh%dm%0.3fs',Hours,Minutes,Seconds);
elseif (Minutes>0)
  msg = sprintf('saffron took %dm%0.3fs',Minutes,Seconds);
else
  msg = sprintf('saffron took %0.3fs',Seconds);
end
logmsg(1,msg);

logmsg(1,'=end=saffron======%s=================\n',datestr(now));

clear global EasySpinLogLevel

return
%=======================================================================


%=======================================================================
% Performs a 2D Inverse Fourier Transformation of a 2D array spc and
% returns the (1:N,1:N) portion of the result. Efficient when N is
% much smaller than the array dimension.
%=======================================================================
% is faster for (128*4, 128*8, 128*16; 256*4, 256*8; 512*4, 512*8)
function td = ifft2dpartial(spc,N,usePartialIFFT)
if (usePartialIFFT)
  % IFFT along one dim, pick 1:N, IFFT along the other, pick 1:N
  td = ifft(spc,[],1);
  td = ifft(td(1:N(1),:),[],2);
  td = td(:,1:N(2));
else
  % IFFT2 the entire array and then pick small block (1:N,1:N)
  td = ifft2(spc);
  td = td(1:N(1),1:N(2));
end