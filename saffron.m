% saffron    Simulate pulse EPR spectra
%
%     [x,S] = saffron(Sys,Exp,Opt)
%     [x,S,out] = saffron(Sys,Exp,Opt)
%
%     [x1,x2,S] = saffron(Sys,Exp,Opt)
%     [x1,x2,S,out] = saffron(Sys,Exp,Opt)
%
%     Sys   ... spin system with electron spin and ESEEM nuclei
%     Exp   ... experimental parameters (time unit us)
%     Opt   ... simulation options
%
%     out:
%       x       ... time or frequency axis (1D experiments)
%       x1, x2  ... time or frequency axis (2D experiments)
%       S       ... simulated signal (ESEEM) or spectrum (ENDOR)
%       out     ... structure with FFT of ESEEM signal

function varargout = saffron(Sys,Exp,Opt)

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

% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

%
% Loop over species and isotopologues
%==================================================================
TwoDim = (isfield(Exp,'Sequence') && ...
  strcmp(Exp.Sequence,'HYSCORE')) || ...
  (isfield(Exp,'Inc') && (max(abs(Exp.Inc))>1));
isENDOR = isfield(Exp,'Sequence') && strcmp(Exp.Sequence,'MimsENDOR');
if ~isfield(Opt,'IsoCutoff'), Opt.IsoCutoff = 1e-4; end

if ~isfield(Sys,'singleiso') || (Sys.singleiso==0)
  
  [SysList,weight] = expandcomponents(Sys,Opt.IsoCutoff);
    
  ysum = 0; % direct domain (TD for ESEEM, FD for ENDOR)
  zsum = 0; % inverse domain (FD for ESEEM)
  for iComponent = 1:numel(SysList)
    if TwoDim
      [x1,x2,y_,out] = saffron(SysList{iComponent},Exp,Opt);
    else
      [x1,y_,out] = saffron(SysList{iComponent},Exp,Opt);
    end
    ysum = ysum + y_*weight(iComponent);
    if ~isENDOR
      zsum = zsum + out.fd*weight(iComponent);
    end
  end
 
  if ~isENDOR
    out.fd = zsum;
    out.td = ysum;
  else
    out.fd = ysum;
  end

  switch nargout
    case 0, % plotting, done below
    case 1, varargout = {ysum};
    case 2, varargout = {x1,ysum};
    case 3,
      if TwoDim
        varargout = {x1,x2,ysum};
      else
        varargout = {x1,ysum,out};
      end
    case 4,
      if TwoDim
        varargout = {x1,x2,ysum,out};
      end
  end

  %===============================================================
  % Plotting
  %===============================================================
  ShowPlots = (nargout==0);
  if ShowPlots
    logmsg(1,'Graphical rendering...');
    clf
    
    if isENDOR
      
      plot(x1,out.fd);
      xlabel('frequency (MHz)');
      ylabel('intensity (arb.units)');
      if isfield(Exp,'tau') && (Exp.tau<1)
        h=title(sprintf('Mims ENDOR, %g mT, tau = %g ns',Exp.Field,1000*Exp.tau));
      else
        h=title(sprintf('Mims ENDOR, %g mT, tau = %g us',Exp.Field,Exp.tau));
      end
      set(h,'FontWeight','b');
      
    else
      
      if ~TwoDim
        
        % Time domain
        subplot(2,1,1);
        PredefinedExperiment = isfield(Exp,'Sequence');
        ExpNames = {'2pESEEM','3pESEEM','4pESEEM','HYSCORE','MimsENDOR'};
        QuadratureSignal = ~PredefinedExperiment;
        if QuadratureSignal
          h = plot(x1,real(out.td),'b',x1,imag(out.td),'r');
          set(h(1),'Color',[0 0 1]);
          set(h(2),'Color',[0.8 0.6 1]);
          legend('Re','Im');
          legend boxoff
        else
          plot(x1,real(out.td));
        end
        axis tight
        xl = xlim;
        xlim([0 xl(2)]);
        yl = ylim;
        ylim(yl+[-1 1]*diff(yl)*0.1);
        
        if PredefinedExperiment
          ExperimentID = strmatch(Exp.Sequence,ExpNames);
          xlb = {'\tau (\mus)','\tau+T (\mus)','T (\mus)','...','frequency (MHz)'};
          xlabel(xlb{ExperimentID});
          ylabel('echo amplitude');
          title([ExpNames{ExperimentID},', TD signal']);
        else
          xlabel('t (\mus)');
          ylabel('echo amplitude (arb.u.)');
          title('User-defined experiment, TD signal');
        end
        set(gca,'Layer','top');
        
        % Frequency domain
        subplot(2,1,2);
        idx = find(out.f==0):length(out.f);
        xf = out.f(idx);
        if QuadratureSignal
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
        if isfield(Sys,'Nucs')
          nuI = larmorfrq(Sys.Nucs,Exp.Field);
          for k=1:numel(nuI)
            line([1 1]*abs(nuI(k)),ylim,'Color',[1 1 1]*0.8);
          end
          h = get(gca,'Children');
        end
        set(gca,'Children',h(end:-1:1));
        
      else
        
        subplot(1,2,1);
        
        pcolor(x1,x2,real(out.td.')); shading flat; axis equal tight;
        set(gca,'Layer','top');
        title('Time domain (real part), baseline corrected');
        xlabel('t_1 (\mus)');
        ylabel('t_2 (\mus)');
        
        subplot(1,2,2);
        
        fx1 = fdaxis(Exp.dt(1),size(out.fd,1));
        if numel(Exp.dt)<2, Exp.dt(2) = Exp.dt(1); end
        fx2 = fdaxis(Exp.dt(2),size(out.fd,2));
        if isfield(Opt,'logplot') && Opt.logplot
          fd = log(out.fd);
          maxfd = max(max(fd));
          fd(fd<maxfd-6) = maxfd-6;
        end
        pcolor(fx1,fx2,out.fd.'); shading flat; axis equal tight
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
      end
      
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
if (Sys.nNuclei==0)
  %error('saffron: There are no nuclear spins in the spin system.');
end

if (Sys.nElectrons>1)
  error('saffron does not support systems with more than one electron spin.');
end
maxNuclei = 40;
if (Sys.nNuclei>maxNuclei)
  error('saffron does not support systems with more than %d nuclei.',maxNuclei);
end

if isfield(Sys,'ExciteWidth')
  error('You gave Sys.ExciteWidth, but it should be Exp.ExciteWidth.');
end

%===================================================================
% Experiment structure
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
if ~isfield(Exp,'T1T2')
  Exp.T1T2 = [0 0];
end
Exp.T1T2(Exp.T1T2==0) = inf;
if numel(Exp.T1T2)~=2
  error('Exp.T1T2 must contain two numbers, T1 and T2 in microseconds.');
end
if any(Exp.T1T2<=0) || any(~isreal(Exp.T1T2))
  error('T1 and T2 in Exp.T1T2 must be positive, in microseconds.');
end

% Pulse sequence
PredefinedExperiment = isfield(Exp,'Sequence');
if PredefinedExperiment
  
  QuadratureSignal = 0;
  
  if isfield(Exp,'Filter')
    error('Exp.Filter can only be used with custom sequences.');
  end

  ExpNames = {'2pESEEM','3pESEEM','4pESEEM','HYSCORE','MimsENDOR'};
  ExperimentID = strmatch(Exp.Sequence,ExpNames);
  if isempty(ExperimentID)
    error('Exp.Sequence ''%s'' not recognized.',Exp.Sequence);
  end
  if numel(ExperimentID)>1, error('Ambiguous sequence name.'); end
  logmsg(1,'Sequence: %s',ExpNames{ExperimentID});
  
  if isfield(Exp,'tp')
    if any(Exp.tp~=0)
      error('You cannot use predefined sequences (Exp.Sequence) with real pulses (Exp.tp).');
    end
  end
    
  isENDOR = false;
  switch ExperimentID
    case 1 % 2pESEEM
      nIntervals = 2; nDimensions = 1; IncSchemeID = 2; nPathways = 1; pulseprefactor = +1/2;
    case 2 % 3pESEEM
      nIntervals = 3; nDimensions = 1; IncSchemeID = 1; nPathways = 2; pulseprefactor = +1/8;
    case 3 % 4pESEEM
      nIntervals = 4; nDimensions = 1; IncSchemeID = 2; nPathways = 2; pulseprefactor = -1/8;
    case 4 % HYSCORE
      nIntervals = 4; nDimensions = 2; IncSchemeID = 11; nPathways = 2; pulseprefactor = -1/8;
    case 5 % Mims ENDOR
      isENDOR = true;
      nIntervals = 3; nDimensions = 1; IncSchemeID = 0; nPathways = 2; pulseprefactor = +1/8;
  end
  
  if ~isfield(Exp,'tau')
    if ExperimentID==1
      Exp.tau = 0;
    else
      error('Exp.tau is missing.');
    end
  end
  if (ExperimentID>1) && all(Exp.tau==0)
    error('Exp.tau must be larger than 0.');
  end
    
else
  
  % User-specified pulse sequence -----------------------------------------
  logmsg(1,'User-specified pulse experiment.');
  ExperimentID = -1;
  isENDOR = 0;
  
  if any(~isinf(Exp.T1T2))
    error('Sorry, T1 and T2 for custom sequences not supported.');
  end
  
  if isfield(Exp,'Phase')
    QuadratureSignal = 1;
  else
    QuadratureSignal = 0;
  end
  
  if ~isfield(Exp,'Flip')
    error('For a pulse experiment, give either Exp.Sequence or Exp.Flip/Exp.Inc.');
  end
  nIntervals = numel(Exp.Flip);
  
  % Incrementation scheme
  if ~isfield(Exp,'Inc')
    error('Exp.Inc with the incrementation scheme is missing.');
  end
  if numel(Exp.Inc)~=nIntervals
    error('Exp.Inc must contain the same number of elements as Exp.Flip.');
  end  
  IncScheme = Exp.Inc(Exp.Inc~=0);
  IncSchemeID = 0;
  switch numel(IncScheme)
    case 1
      if (IncScheme==1), IncSchemeID = 1; end
    case 2
      if all(IncScheme==[1 1]), IncSchemeID = 2;
      elseif all(IncScheme==[1 -1]), IncSchemeID = 3;
      elseif all(IncScheme==[1 2]), IncSchemeID = 11;
      end
    case 3
      if all(IncScheme==[1 2 1]), IncSchemeID = 12;
      elseif all(IncScheme==[1 2 2]), IncSchemeID = 13;
      elseif all(IncScheme==[1 1 2]), IncSchemeID = 14;
      end
    case 4
      if all(IncScheme==[1 2 2 1]), IncSchemeID = 15;
      elseif all(IncScheme==[1 2 -2 1]), IncSchemeID = 16;
      elseif all(IncScheme==[1 1 2 2]), IncSchemeID = 17;
      end
  end
  if (IncSchemeID==0)
    error('Unrecognized incrementation scheme.');
  end
  nDimensions = max(abs(Exp.Inc));
  
  if ~isfield(Exp,'t')
    Exp.t = zeros(1,nIntervals);
  end
  if numel(Exp.t)~=nIntervals
    error('Exp.t must contain the same number of elements as Exp.Flip');
  end
  if all(Exp.t==0) && (numel(IncScheme)<nIntervals)
    logmsg(0,'Some delays are zero, but are not incremented!');
  end

end

% dt
if ~isENDOR
  if ~isfield(Exp,'dt')
    error('Exp.dt is missing.');
  end
  if numel(Exp.dt)==1
    Exp.dt = Exp.dt*ones(1,nDimensions);
  elseif numel(Exp.dt)~=nDimensions
    error('Exp.dt needs either 1 or %d elements, one per dimension. You gave %d.',nDimensions,numel(Exp.dt));
  end
end

% nPoints
if ~isfield(Exp,'nPoints')
  if isENDOR
    Exp.nPoints = 1001;
  else
    if (nDimensions==1)
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

% tp
% Determine whether to use ideal pulse theory
if ~isfield(Exp,'tp')
  Exp.tp = zeros(1,nIntervals);
else
  if numel(Exp.tp)~=nIntervals
    error('Exp.tp contains a wrong number of elements.');
  end
end
idealPulse = (Exp.tp==0);
realPulse = ~idealPulse;

% excitation width
if any(realPulse)
  if isfield(Exp,'ExciteWidth')
    %error('Cannot Exp.ExciteWidth and real pulses (Exp.tp) at the same time.')
  end
end


% User-defined experiment: pathways etc
%--------------------------------------------------------------------------

if ~PredefinedExperiment

  % Validate pulse sequence
  Valid = 1;
  for iInterval = 1:nIntervals
    t_start = Exp.t(iInterval);
    Valid = (t_start>=0);
    iDim = abs(Exp.Inc(iInterval));
    if (iDim>0)
      t_end = t_start + ...
        sign(Exp.Inc(iInterval))*(Exp.nPoints(iDim)-1)*Exp.dt(iDim);
      if (t_end<0), Valid = 0; end
    end
    if ~Valid, break; end
  end
  if (~Valid)
    error('Pulse sequence gives negative length after pulse %d.',iInterval);
  end

  % Determine pathways contributing to the echo
  logmsg(1,'  determining pathways contributing to the echo');
  if isfield(Exp,'Pathways')
    code0('ab+-') = [1 2 3 4];
    pathwayList = code0(Exp.Pathways);
    if size(pathwayList,2)~=nIntervals
      error('Exp.Pathways contains pathways longer or shorter than the pulse sequence.');
    end
  else
    pathwayList = sf_pathways(Exp);
  end
  
  nPathways = size(pathwayList,1);
  if (nPathways==0)
    error('Sorry, no focused echo with this sequence and timings.');
  end
  
  % Apply coherence filters ---------------------
  if ~isfield(Exp,'Filter')
    Exp.Filter = [];
  end
  if ~isempty(Exp.Filter)
    logmsg(1,'  applying user-supplied coherence filters');
    if ~ischar(Exp.Filter)
      error('Exp.Filter must be a string containing ''0'', ''1'', ''a'', ''b'', ''+'', ''-'' and/or ''.''.');
    end
    if (numel(Exp.Filter)~=nIntervals)
      error('Exp.Filter must have %d elements instead of the given %d elements.',nIntervals,numel(Exp.Filter));
    end
    for iInt = 1:nIntervals
      switch Exp.Filter(iInt)
        case '0', keep = (pathwayList(:,iInt)==1) | (pathwayList(:,iInt)==2);
        case '1', keep = (pathwayList(:,iInt)==3) | (pathwayList(:,iInt)==4);
        case 'a', keep = (pathwayList(:,iInt)==1);
        case 'b', keep = (pathwayList(:,iInt)==2);
        case '+', keep = (pathwayList(:,iInt)==3);
        case '-', keep = (pathwayList(:,iInt)==4);
        case '.', keep = (pathwayList(:,iInt)~=0);
        case '*', keep = (pathwayList(:,iInt)~=0);
        otherwise
          error('Exp.Filter(%d) = ''%s'' is invalid. Use ''0'', ''1'', ''a'', ''b'', ''+'', ''-'' or ''.''.',iInt,Exp.Filter(iInt));
      end
      pathwayList = pathwayList(keep,:);
    end
    nPathways = size(pathwayList,1);
    if (nPathways==0)
      error('Exp.Filter is too restrictive: no echo at detection point left after applying the filter.');
    end
  else
    logmsg(1,'  no user-supplied coherence filters');
  end

  [idxFreeL,idxFreeR,idxPulseL,idxPulseR] = pathwayparser(pathwayList);

  % Phase is in multiples of pi/2
  if ~isfield(Exp,'Phase')
    Exp.Phase = ones(1,nIntervals);
  end

  % compute all ideal pulse transfer factors
  pulseprefactor = ones(1,nPathways);
  for iPulse = 1:numel(Exp.t)
    if idealPulse(iPulse)
      theta = Exp.Flip(iPulse)*pi/2;
      c = cos(theta/2);
      s = sin(theta/2);
      for iPathway = 1:nPathways
        switch idxPulseL(iPathway,iPulse)
          case 1, pL = c;
          case 2, pL = c;
          %case 3, pL = -1i*s*exp(-1i*Exp.Phase(iPulse));
          %case 4, pL = -1i*s*exp(+1i*Exp.Phase(iPulse));
          case 3, pL = -1i*s*(-1i)^Exp.Phase(iPulse);
          case 4, pL = -1i*s*(+1i)^Exp.Phase(iPulse);
        end
        switch idxPulseR(iPathway,iPulse)
          case 1, pR = c;
          case 2, pR = c;
          %case 3, pR = +1i*s*exp(+1i*Exp.Phase(iPulse));
          %case 4, pR = +1i*s*exp(-1i*Exp.Phase(iPulse));
          case 3, pR = +1i*s*(+1i)^Exp.Phase(iPulse);
          case 4, pR = +1i*s*(-1i)^Exp.Phase(iPulse);
        end
        pulseprefactor(iPathway) = pL*pulseprefactor(iPathway)*pR;
      end
    end
  end
     
  % Remove pathways with pulseprefactor zero
  rmv = abs(pulseprefactor)<1e-6;
  pathwayList(rmv,:) = [];
  idxFreeL(rmv,:) = [];
  idxFreeR(rmv,:) = [];
  idxPulseL(rmv,:) = [];
  idxPulseR(rmv,:) = [];
  pulseprefactor(rmv) = [];

  nPathways = size(pathwayList,1);

  idxIncL = idxFreeL(:,Exp.Inc~=0);
  idxIncR = idxFreeR(:,Exp.Inc~=0);

  if (EasySpinLogLevel>0)
    logmsg(1,'  Pathways and prefactors');
    Str = 'ab+-';
    for iPathway = 1:nPathways
      logmsg(1,'  %s  (%+4.3f, %+4.3f)',Str(pathwayList(iPathway,:)),real(pulseprefactor(iPathway)),imag(pulseprefactor(iPathway)));
    end
  end

end




if isENDOR
  
  if all(Sys.lwEndor==0)
    error('Positive ENDOR line width in Sys.lwEndor needed.');
  end
  Sys.lwEndor = Sys.lwEndor(1);

  if ~isfield(Exp,'Range')
    error('Frequency range (Exp.Range) must be given for an ENDOR experiment.');
  end

end


% Relaxation time constants
if ~isfield(Exp,'T1'), Exp.T1 = 0; end
if ~isfield(Exp,'T2'), Exp.T2 = 0; end


if PredefinedExperiment
  switch ExperimentID
    case 2 % 3pESEEM
      if ~isfield(Exp,'T'), Exp.T = 0; end
    case 5 % Mims ENDOR
      if ~isfield(Exp,'T'), Exp.T = 0; end
    case 3 % 4pESEEM
      if ~isfield(Exp,'T'), Exp.T = 0; end
    case 4 % HYSCORE
      if ~isfield(Exp,'t1'), Exp.t1 = 0; end
      if ~isfield(Exp,'t2'), Exp.t2 = 0; end
      if numel(Exp.t1)>1
        error('Exp.t1 must a single positive number, in units of us.');
      end
      if numel(Exp.t2)>1
        error('Exp.t2 must a single positive number, in units of us.');
      end
      if (Exp.t1~=Exp.t2)
        fprintf('Exp.t1 and Exp.t2 are not identical, so the resulting spectrum might be asymmetric.\n');
      end
  end
end


% Pick nuclei to be included in the computation
shfNuclei = 1:Sys.nNuclei;
if all(idealPulse) && isfield(Exp,'ExciteWidth')
  idxStrongNuclei = max(abs(Sys.A),[],2)>Exp.ExciteWidth;
  shfNuclei(idxStrongNuclei) = [];
end

TwoElectronManifolds = (Sys.nElectrons==1) && (Sys.S==1/2) && ...
  (numel(shfNuclei)==Sys.nNuclei);

OrientationSelection = isfield(Exp,'mwFreq');
if (OrientationSelection)
  logmsg(1,'Microwave frequency given: orientation selection is on.');
  if ~isfield(Exp,'ExciteWidth')
    error('Orientation selection: Exp.ExciteWidth (in MHz) missing. It should be about the inverse of the first pulse length (100MHz for 10ns). If you don''t want orientation selection, set it to a very large number (1e6) or remove the microwave frequency.');
  end
else
  logmsg(1,'no orientation selection (infinite bandwidth).');
end

if isfield(Exp,'ExciteWidth')
  if ~isfield(Exp,'mwFreq')
    error('Exp.ExciteWidth is given, but Exp.mwFreq is missing. Please give Exp.mwFreq.');
  end
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
if ~isfield(Opt,'Transitions')
  Opt.Transitions = [];
end

% Nuclei: which nuclei to include in the simulation
if isfield(Opt,'Nuclei')
  if isempty(Opt.Nuclei)
    error('Opt.Nuclei must contain the indices of nuclei to include in the simulation.');
  end
  if any(Opt.Nuclei>Sys.nNuclei) || any(Opt.Nuclei<0)
    error('Index in Opt.Nuclei out of range.');
  end
  shfNuclei = Opt.Nuclei;
else
  % Pick nuclei to be included in the ESEEM/ENDOR computation
  shfNuclei = 1:Sys.nNuclei;
  if all(idealPulse) && isfield(Exp,'ExciteWidth')
    idxStrongNuclei = max(abs(Sys.A),[],2)>Exp.ExciteWidth;
    shfNuclei(idxStrongNuclei) = [];
  end
end

% Expansion factor: determines size of spectral buffer
if ~isfield(Opt,'Expand')
  if (nDimensions==1), Opt.Expand = 4; else Opt.Expand = 2; end
end
maxExpand = 8;
if (numel(Opt.Expand)~=1) || (Opt.Expand<0) || (Opt.Expand>maxExpand) || rem(Opt.Expand,1)
  error('Opt.Expand must be an integer between 0 and %d.',maxExpand);
end

% Number of knots: determines number of orientations
if ~isfield(Opt,'nKnots'), Opt.nKnots = 30+1; end
if numel(Opt.nKnots)>1
  error('Only one number allowed in Opt.nKnots. saffron does not support interpolation.');
end
if (Opt.nKnots<7)
  error('Opt.nKnots must be at least 7. You gave %d.',Opt.nKnots);
end

% ProductRule: determines whether product rule is used or not
if ~isfield(Opt,'ProductRule'), Opt.ProductRule = 0; end
if (Sys.nNuclei==1), Opt.ProductRule = 0; end
%if (isENDOR), Opt.ProductRule = 1; end

if (any(realPulse) && Opt.ProductRule)
  error('saffron: Cannot apply product rule and real pulses at the same time.');
end

if isfield(Opt,'Output'), logmsg(0,'saffron does not support Opt.Output.'); end
SeparateTransitions = false;

if ~isfield(Opt,'OriThreshold'), Opt.OriThreshold = 0.005; end

if ~isfield(Opt,'Window'),
  if (nDimensions==1), Opt.Window = 'ham+'; else Opt.Window = 'ham+'; end
end

if ~isfield(Opt,'ZeroFillFactor');
  Opt.ZeroFillFactor = 2;
end

%DataProcessing = ShowPlots;
DataProcessing = 1;

% undocumented options
if ~isfield(Opt,'logplot'), Opt.logplot = 0; end
if ~isfield(Opt,'PartialIFFT'), Opt.PartialIFFT = 1; end
if ~isfield(Opt,'TimeDomain'), Opt.TimeDomain = 0; end


%==========================================================================
% Symmetry determination and orientational grid.
%==========================================================================
p_symandgrid;

p_crystalorientations;

logmsg(1,'-Hamiltonians------------------------------------------');

%=====================================================================
% Compute electronic Hamiltonian
%=====================================================================

logmsg(1,'setting up electronic Hamiltonian...');
if (TwoElectronManifolds)
  % not needed, S=1/2 electronic Hamiltonian can be solved analytically
else
  coreSys = nucspinrmv(Sys,shfNuclei);
  coreSys.processed = 0;
  coreSys.lw = 0;
  % Operators for constructing Hamiltonian
  [F,Gx,Gy,Gz] = sham(coreSys);
  % Operators for computing <i|S|i>
  Sx = sop(coreSys,1,1); % works only for one electron spin
  Sy = sop(coreSys,1,2);
  Sz = sop(coreSys,1,3);
end
%=====================================================================


%=====================================================================
% Compute nuclear spin Hamiltonians
%=====================================================================
logmsg(1,'computing nuclear spin sub-Hamiltonians...');

SeparateSpaces = Opt.ProductRule;
if SeparateSpaces
  logmsg(1,'  separate subspace for each nucleus');
  for iNuc = shfNuclei % only those explicitly included
    [Ix,Iy,Iz] = sop(Sys.I(iNuc),'x','y','z');
    
    % Nuclear Zeeman -------------------------------------------------
    gamm = -Sys.gn(iNuc)*nmagn/1e3/planck/1e6;
    Space(iNuc).Hnz1 = gamm*Ix;
    Space(iNuc).Hnz2 = gamm*Iy;
    Space(iNuc).Hnz3 = gamm*Iz;
    
    % Hyperfine ------------------------------------------------------
    if Sys.fullA
      A = Sys.A((iNuc-1)*3+(1:3),:);
    else
      Ra = eye(3);
      if isfield(Sys,'Apa'), Ra = erot(Sys.Apa(iNuc,:)); end
      A = Ra*diag(Sys.A(iNuc,:))*Ra';
    end
    Space(iNuc).Hhf1 = A(1,1)*Ix + A(1,2)*Iy + A(1,3)*Iz;
    Space(iNuc).Hhf2 = A(2,1)*Ix + A(2,2)*Iy + A(2,3)*Iz;
    Space(iNuc).Hhf3 = A(3,1)*Ix + A(3,2)*Iy + A(3,3)*Iz;
    
    % Nuclear quadrupole ---------------------------------------------
    if (Sys.I(iNuc)>=1)
      Rq = erot(Sys.Qpa(iNuc,:));
      Q = Rq*diag(Sys.Q(iNuc,:))*Rq'; Q = (Q + Q.')/2;
      Space(iNuc).Hnq = ...
        Ix*(Q(1,1)*Ix + Q(1,2)*Iy + Q(1,3)*Iz) + ...
        Iy*(Q(2,1)*Ix + Q(2,2)*Iy + Q(2,3)*Iz) + ...
        Iz*(Q(3,1)*Ix + Q(3,2)*Iy + Q(3,3)*Iz);
    else
      Space(iNuc).Hnq = 0;
    end
  end
    
else
  logmsg(1,'  complete nuclear state space');
  
  Space.Hnq = 0;
  Space.Hhf1 = 0; Space.Hhf2 = 0; Space.Hhf3 = 0;
  Space.Hnz1 = 0; Space.Hnz2 = 0; Space.Hnz3 = 0;
  I = Sys.I(shfNuclei);
  for ishfNuc = 1:numel(shfNuclei)
    Ix = sop(I,ishfNuc,1);
    Iy = sop(I,ishfNuc,2);
    Iz = sop(I,ishfNuc,3);
    iNuc = shfNuclei(ishfNuc);
    
    % Nuclear Zeeman -------------------------------------------------
    gamm = -Sys.gn(iNuc)*nmagn/1e3/planck/1e6;
    Space.Hnz1 = Space.Hnz1 + gamm*Ix;
    Space.Hnz2 = Space.Hnz2 + gamm*Iy;
    Space.Hnz3 = Space.Hnz3 + gamm*Iz;
    
    % Hyperfine ------------------------------------------------------
    if Sys.fullA
      A = Sys.A((iNuc-1)*3+(1:3),:);
    else
      Ra = eye(3);
      if isfield(Sys,'Apa'), Ra = erot(Sys.Apa(iNuc,:)); end
      A = Ra*diag(Sys.A(iNuc,:))*Ra';
    end
    Space.Hhf1 = Space.Hhf1 + A(1,1)*Ix + A(1,2)*Iy + A(1,3)*Iz;
    Space.Hhf2 = Space.Hhf2 + A(2,1)*Ix + A(2,2)*Iy + A(2,3)*Iz;
    Space.Hhf3 = Space.Hhf3 + A(3,1)*Ix + A(3,2)*Iy + A(3,3)*Iz;
    
    % Nuclear quadrupole ---------------------------------------------
    if (Sys.I(iNuc)>=1)
      if Sys.fullQ
        Q = Sys.Q((iNuc-1)*3+(1:3),:);
      else
        Rq = eye(3);
        if isfield(Sys,'Qpa'), Rq = erot(Sys.Qpa(iNuc,:)); end
        Q = Rq*diag(Sys.Q(iNuc,:))*Rq';
        Q = (Q + Q.')/2;
      end
      Hnq_ = ...
        Ix*(Q(1,1)*Ix + Q(1,2)*Iy + Q(1,3)*Iz) + ...
        Iy*(Q(2,1)*Ix + Q(2,2)*Iy + Q(2,3)*Iz) + ...
        Iz*(Q(3,1)*Ix + Q(3,2)*Iy + Q(3,3)*Iz);
    else
      Hnq_ = 0;
    end
    Space.Hnq = Space.Hnq + Hnq_;
  end
  
end
nSubSpaces = numel(Space);
logmsg(1,'  %d nuclei, %d subspaces',numel(shfNuclei),nSubSpaces);
%=====================================================================



%=====================================================================
% Orientation pre-selection factors
%=====================================================================
OrientationPreSelection = OrientationSelection && TwoElectronManifolds;
if OrientationPreSelection
  logmsg(1,'pre-computing orientation selecton from g tensor alone...');
  logmsg(1,'  S=1/2: using simple g tensor/HStrain model');
  logmsg(1,'  excitation width (MHz): %g',Exp.ExciteWidth);
  logmsg(1,'  HStrain (MHz): %g %g %g',Sys.HStrain(1),Sys.HStrain(2),Sys.HStrain(3));
  % g values for all orientations
  zLab = ang2vec(Orientations(:,1),Orientations(:,2));
  geff = diag(Sys.g)*zLab;
  geff = sqrt(sum(geff.^2,1));
  % resonance frequencies for all orientations
  nu = geff*bmagn*Exp.Field/1e3/planck/1e6; % MHz
  % line widths for all orientations
  lw = diag(Sys.HStrain)*zLab;
  lw2 = (sum(lw.^2,1)) + Exp.ExciteWidth^2;
  gOriSelWeight = exp(-(nu-1000*Exp.mwFreq).^2./lw2);
end

%=====================================================================
% Preparation
%=====================================================================
logmsg(1,'Preparation...');

if isENDOR
  
  rf = linspace(Exp.Range(1),Exp.Range(2),Exp.nPoints);
  Template.x0 = 5e4;
  Template.lw = Template.x0/2.5; %<1e-8 at borders for Harmonic = -1
  Template.y = gaussian(0:2*Template.x0-1,Template.x0,Template.lw,-1);
  Template.y = Template.y*(rf(2)-rf(1))/Sys.lwEndor;
  %plot(Template.y);
  
  endorspc = zeros(1,Exp.nPoints);
  endoroffset = 0;

else

  if Opt.TimeDomain
    logmsg(1,'  time domain simulation');
  else
    logmsg(1,'  frequency domain simulation');
    % Prepare for binning method
    ExpansionFactor = 2.^Opt.Expand;
    nPointsF = ExpansionFactor*Exp.nPoints;
    if numel(nPointsF)==2
      logmsg(1,'  %dx%d points, expand x%d -> %dx%d points',...
        Exp.nPoints(1),Exp.nPoints(2),ExpansionFactor,nPointsF(1),nPointsF(2));
    else
      logmsg(1,'  %d points, expand x%d -> %d points',Exp.nPoints,ExpansionFactor,nPointsF);
    end
    % allocate array(s)
    if (nDimensions==1), siz = [1, nPointsF]; else siz = nPointsF; end
    if Opt.ProductRule
      for iP = 1:nPathways
        for iS = 1:nSubSpaces
          pathwaybuff{iP,iS} = zeros(siz); 
          pathwaybuff{iP,iS}(1) = 1e-300i; % make sure it's complex
        end
      end
    else
      buff = zeros(siz);
      buff(1) = 1e-300i; % make sure it's complex
    end
  end
  
  totaltd = 0;

end
%=====================================================================





%=====================================================================
% Orientation loop
%=====================================================================
logmsg(1,'Looping over %d orientations...',nOrientations);

% Prepare offsets
if any(realPulse)
  if ~isfield(Opt,'nOffsets');
    Opt.nOffsets = 291;
  end
  if ~isfield(Opt,'lwOffset')
    Opt.lwOffset = 100;
  end
  offsets = linspace(-1,1,Opt.nOffsets)*Opt.lwOffset*2;
  offsetWeight = exp(-(offsets/Opt.lwOffset).^2);
  offsetWeight = offsetWeight/sum(offsetWeight);
else
  offsets = 0;
  offsetWeight = 1;
  Opt.nOffsets = 1;
end

if (TwoElectronManifolds)
  g = Sys.g;
  if ~Sys.fullg
    R = erot(Sys.gpa);
    g = R*diag(g)*R';
  end
end

nSkippedOrientations = 0;
for iOri = 1:nOrientations
  
  
  % Get magnetic field orientation
  %------------------------------------------------------------------
  [xLab_M,yLab_M,zLab_M] = erot(Orientations(iOri,:));
  % xLab_M, yLab_M, zLab_M represented in the molecular frame

  % Compute electronic Hamiltonian, energies and <S>
  %------------------------------------------------------------------
  if (TwoElectronManifolds)    
    
    if (OrientationSelection)
      OriSelWeight = gOriSelWeight(iOri);
    else
      OriSelWeight = 1;
    end
    if (OriSelWeight<Opt.OriThreshold)
      nSkippedOrientations = nSkippedOrientations + 1;
      continue
    end
    
    quantizationAxis = g.'*zLab_M(:);
    quantizationAxis = quantizationAxis/norm(quantizationAxis);
    Manifold(1).S = -0.5*quantizationAxis;
    Manifold(2).S = +0.5*quantizationAxis;

    Transitions = [1 2];
    nTransitions = 1;
    ManifoldsInvolved = [1 2];

  else

    % automatic transition selection
    %------------------------------------------------------------
    if isempty(Opt.Transitions)
      H = F + Exp.Field*(zLab_M(1)*Gx + zLab_M(2)*Gy + zLab_M(3)*Gz);
      [eV,eE] = eig(H);
      SyLab = yLab_M(1)*Sx + yLab_M(2)*Sy + yLab_M(3)*Sz;
      SyLab = abs(eV'*SyLab*eV);
      maxSyy = max(SyLab(:));
      
      if (OrientationSelection)
        eE = real(diag(eE));
        dE = eE(:,ones(1,length(eE)));
        dE = (dE-dE') - Exp.mwFreq*1e3; % MHz
        excitationAmplitude = exp(-(dE/Exp.ExciteWidth).^2);
        SyLab = SyLab.*excitationAmplitude;
      end
      
      SyLab(SyLab<Opt.OriThreshold*maxSyy) = 0;
      [v,u,OriSelWeight] = find(tril(SyLab,-1));
    else
      u = Opt.Transitions(:,1);
      v = Opt.Transitions(:,2);
      OriSelWeight = ones(1,length(u));
    end
    Transitions = [u,v];
    nTransitions = size(Transitions,1);

    if (nTransitions==0)
      nSkippedOrientations = nSkippedOrientations + 1;
      continue
    end

    % computation of <S> for all manifolds involved
    %------------------------------------------------------------
    ManifoldsInvolved = zeros(1,length(H));
    ManifoldsInvolved(u) = 1;
    ManifoldsInvolved(v) = 1;
    ManifoldsInvolved = find(ManifoldsInvolved);
    for iM = ManifoldsInvolved
      vec = eV(:,iM);
      Manifold(iM).S = real([vec'*Sx*vec; vec'*Sy*vec; vec'*Sz*vec]);
    end
  end
    
  logmsg(2,'orientation %d of %d: %d transitions',iOri,nOrientations,nTransitions);
  
  % Compute and diagonalize nuclear Hamiltonians
  %----------------------------------------------------------------------
  for iSpace = 1:nSubSpaces
    Hnuc = Exp.Field*(zLab_M(1)*Space(iSpace).Hnz1 + zLab_M(2)*Space(iSpace).Hnz2 + zLab_M(3)*Space(iSpace).Hnz3) + Space(iSpace).Hnq;
    for iM = ManifoldsInvolved
      S = Manifold(iM).S;
      H = Hnuc + S(1)*Space(iSpace).Hhf1 + ...
                 S(2)*Space(iSpace).Hhf2 + ...
                 S(3)*Space(iSpace).Hhf3;
      [VV,EE] = eig((H+H')/2);
      Manifold(iM).V{iSpace} = VV;
      Manifold(iM).E{iSpace} = real(diag(EE));
    end
  end
  
  
  % Loop over all excited EPR transitions
  %----------------------------------------------------------------------
  for iT = 1:nTransitions
    b = Transitions(iT,1); % lower manifold
    a = Transitions(iT,2); % upper manifold

    % Loop over all subspaces
    for iSpace = 1:nSubSpaces
      Ea = Manifold(a).E{iSpace};
      Eb = Manifold(b).E{iSpace};
      Ma = Manifold(a).V{iSpace};
      Mb = Manifold(b).V{iSpace};
      M = Ma'*Mb;       % <a|b> overlap matrix
      Mt = M';
      nStates = length(Ea);
      eyeN = eye(nStates);
      if any(realPulse)
        idxa = 1:nStates;
        idxb = idxa + nStates;
      end
      
      for iOffset = 1:Opt.nOffsets

        % Pulse propagators
        for iInt = 1:nIntervals
          if realPulse(iInt)
            nu1 = (pi/2*Exp.Flip(iInt))/Exp.tp(iInt)/2/pi;
            Hpulse = [diag(Ea+offsets(iOffset)/2), +M*nu1/2i; ...
              -Mt*nu1/2i diag(Eb-offsets(iOffset)/2)];
            FullPuls = expm(-2i*pi*Exp.tp(iInt)*Hpulse);
            P{iInt,1} = FullPuls(idxa,idxa);
            P{iInt,2} = FullPuls(idxb,idxb);
            P{iInt,3} = FullPuls(idxa,idxb);
            P{iInt,4} = FullPuls(idxb,idxa);
          end
        end

        if isENDOR
          prefactor = Weights(iOri)*OriSelWeight(iT);
        else
          if ~Opt.ProductRule
            prefactor = Weights(iOri)*OriSelWeight(iT);
          else
            prefactor = 1;
          end
          if ~Opt.TimeDomain
            if Opt.ProductRule
              for iP = 1:nPathways
                pathwaybuff{iP,iSpace} = zeros(siz);
                pathwaybuff{iP,iSpace}(1) = 1e-300i; % make sure it's complex
              end
            end
          end
        end
        prefactor = prefactor*offsetWeight(iOffset);

        % Compute peaks and amplitudes
        if ~PredefinedExperiment

          % General method ------------------------------------------
          increments = Exp.Inc;
          BlockL = cell(0);
          BlockR = cell(0);
          for iPathway = 1:nPathways
            iBlock = 0;
            Left = prefactor*pulseprefactor(iPathway);
            Right = 1;
            for iInt = 1:nIntervals
              % Pulse propagator
              if idealPulse(iInt)
                switch idxPulseL(iPathway,iInt)
                  case 3, Left = M*Left;
                  case 4, Left = Mt*Left;
                end
                switch idxPulseR(iPathway,iInt)
                  case 3, Right = Right*Mt;
                  case 4, Right = Right*M;
                end
              else
                Left = P{iInt,idxPulseL(iPathway,iInt)}*Left;
                Right = Right*P{iInt,idxPulseR(iPathway,iInt)}';
              end
              
              % Free evolution propagator
              t = Exp.t(iInt);
              if (t>0)
                if (idxFreeL(iPathway,iInt)==1), E = Ea; else E = Eb; end
                Left = diag(exp(-2i*pi*E*t)) * Left;
                if (idxFreeR(iPathway,iInt)==1), E = Ea; else E = Eb; end
                Right = Right * diag(exp(+2i*pi*E*t));
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
                    BlockL{iBlock} = eyeN;
                  else
                    BlockL{iBlock} = Left;
                  end
                  if numel(Right)==1
                    BlockR{iBlock} = eyeN;
                  else
                    BlockR{iBlock} = Right;
                  end
                end
                iBlock = iBlock + 1;
                Left = 1;
                Right = 1;
              end
              
            end % iInt

            if increments(end)~=0, D = M; else D = Right*M*Left; end

            % Acumulate peaks / generated time domain
            if Opt.TimeDomain
              if Opt.ProductRule
                error('Product rule for user-defined experiments not implemented.');
              else
                totaltd = totaltd + ...
                  sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,...
                  idxIncL(iPathway,:),idxIncR(iPathway,:),...
                  Ea,Eb,G,D,BlockL{:},BlockR{:});
              end
            else
              if Opt.ProductRule
                error('Product rule for user-defined experiments not implemented.');
              else
                sf_peaks(IncSchemeID,buff,Exp.dt,...
                  idxIncL(iPathway,:),idxIncR(iPathway,:),...
                  Ea,Eb,G,D,BlockL{:},BlockR{:});
              end
            end

          end % iPathway loop

        else

          switch (ExperimentID)

            % Mims ENDOR ------------------------------------------------
            case 5

              Q = exp(-2i*pi*Ea*Exp.tau)*exp(-2i*pi*Eb*Exp.tau)';
              PG = prefactor*Q.*M; PD = conj(Q).*M;
              G1 = PG*Mt; D1 = PD*Mt; G2 = Mt*PG; D2 = Mt*PD;

              if (Exp.T~=0)
                S1 = (exp(-2i*pi*Ea*Exp.T)*exp(-2i*pi*Ea*Exp.T)').*S1;
                S2 = (exp(-2i*pi*Eb*Exp.T)*exp(-2i*pi*Eb*Exp.T)').*S2;
              end

              % echo signals in absence of RF pulse (off resonant)
              off1 = trace(G1*D1);
              off2 = trace(G2*D2);
              % loop only over nuclear sublevel pairs with Delta mI = 1
              for j=1:nStates-1
                for i=j+1
                  % > RF pulse approximation: apply pi pulse two two-level subsystem
                  % R = [0 -1; 1 0];
                  %S1_ = S1; S1_([j i],[j i]) = R*S1_([j i],[j i])*R';
                  %S2_ = S2; S2_([j i],[j i]) = R*S2_([j i],[j i])*R';
                  % > RF pulse approximation: only swap diagonal elements ii and jj
                  G1_ = G1; q = G1_(i,i); G1_(i,i) = G1_(j,j); G1_(j,j) = q;
                  G2_ = G2; q = G2_(i,i); G2_(i,i) = G2_(j,j); G2_(j,j) = q;
                  ampl = [off1-trace(G1_*D1), off2-trace(G2_*D2)];
                  freq = [Ea(i)-Ea(j), Eb(i)-Eb(j)];
                  %idx = fix(Exp.nPoints*(freq-rf(1))/(rf(end)-rf(1)))+1;
                  %endorspc(idx) = endorspc(idx) + ampl;
                  %endorspc = endorspc + ...
                  %  lisum1i(Template.y,Template.x0,Template.lw,freq,ampl,Sys.lwEndor*[1 1],rf);
                  if SeparateTransitions
                    endorspc(iT,:) = ampl(1)*exp(-((rf-freq(1))/Sys.lwEndor).^2);
                    endorspc(iT,:) = ampl(2)*exp(-((rf-freq(2))/Sys.lwEndor).^2);
                  else
                    endorspc = endorspc + ampl(1)*exp(-((rf-freq(1))/Sys.lwEndor).^2);
                    endorspc = endorspc + ampl(2)*exp(-((rf-freq(2))/Sys.lwEndor).^2);
                  end
                  endoroffset = endoroffset + off1 + off2;
                end
              end

            % HYSCORE ---------------------------------------------------
            case 4

              if ~all(idealPulse)
                error('Pre-defined HYSCORE with real pulses not supported.');
              end
              
              Q = exp(-2i*pi*Ea*Exp.tau)*exp(-2i*pi*Eb*Exp.tau)';
              PG = prefactor*Q.*M;
              PD = conj(Q).*M;
              G1 = PG*Mt; G2 = Mt*PG;
              D1 = Mt*PD; D2 = PD*Mt;
              if (Exp.t1~=0)
                q = exp(-2i*pi*Exp.t1*Ea); G1 = (q*q').*G1;
                q = exp(-2i*pi*Exp.t1*Eb); G2 = (q*q').*G2;
              end
              if (Exp.t2~=0)
                q = exp(+2i*pi*Exp.t2*Eb); D1 = (q*q').*D1;
                q = exp(+2i*pi*Exp.t2*Ea); D2 = (q*q').*D2;
              end

              if Opt.TimeDomain
                if Opt.ProductRule
                  pathwaytd{1,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Mt,M);
                  pathwaytd{2,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,Mt,M);
                else
                  totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Mt,M);
                  totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,Mt,M);
                end
              else
                if Opt.ProductRule
                  sf_peaks(IncSchemeID,pathwaybuff{1,iSpace},Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Mt,M);
                  sf_peaks(IncSchemeID,pathwaybuff{2,iSpace},Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,M,Mt);
                else
                  sf_peaks(IncSchemeID,buff,Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Mt,M);
                  sf_peaks(IncSchemeID,buff,Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,M,Mt);
                end
              end

            % 4pESEEM ---------------------------------------------------
            case 3

              if ~all(idealPulse)
                error('Pre-defined 4p-ESEEM with real pulses not supported.');
              end
              
              Q = exp(-2i*pi*Ea*Exp.tau)*exp(-2i*pi*Eb*Exp.tau)';
              PG = prefactor*Q.*M;
              PD = conj(Q).*M;
              G1 = PG*Mt; G2 = Mt*PG;
              D1 = Mt*PD; D2 = PD*Mt;
              if (Exp.T~=0)
                q = exp(-2i*pi*Exp.T*Ea); G1 = (q*q').*G1;
                q = exp(-2i*pi*Exp.T*Eb); G2 = (q*q').*G2;
                q = exp(+2i*pi*Exp.T*Eb); D1 = (q*q').*D1;
                q = exp(+2i*pi*Exp.T*Ea); D2 = (q*q').*D2;
              end

              if Opt.TimeDomain
                if Opt.ProductRule
                  pathwaytd{1,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Mt,M);
                  pathwaytd{2,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,Mt,M);
                else
                  totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Mt,M);
                  totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,Mt,M);
                end
              else
                if Opt.ProductRule
                  sf_peaks(IncSchemeID,pathwaybuff{1,iSpace},Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Mt,M);
                  sf_peaks(IncSchemeID,pathwaybuff{2,iSpace},Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,M,Mt);
                else
                  sf_peaks(IncSchemeID,buff,Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Mt,M);
                  sf_peaks(IncSchemeID,buff,Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,M,Mt);
                end
              end

            % 3pESEEM ------------------------------------------------------------
            case 2

              if ~all(idealPulse)
                error('Pre-defined 3pESEEM not supported for real pulses.');
              end
              
              Q = exp(-2i*pi*Ea*Exp.tau)*exp(-2i*pi*Eb*Exp.tau)';
              PG = prefactor*Q.*M; PD = conj(Q).*M;
              G1 = PG*Mt; D1 = PD*Mt; G2 = Mt*PG; D2 = Mt*PD;
              if (Exp.T~=0)
                q = exp(-2i*pi*Exp.T*Ea); G1 = (q*q').*G1;
                q = exp(-2i*pi*Exp.T*Eb); G2 = (q*q').*G2;
              end

              if Opt.TimeDomain
                if Opt.ProductRule
                  pathwaytd{1,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,1,1,Ea,Eb,G1,D1);
                  pathwaytd{2,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,2,2,Ea,Eb,G2,D2);
                else
                  totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,1,1,Ea,Eb,G1,D1);
                  totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,2,2,Ea,Eb,G2,D2);
                end
              else
                if Opt.ProductRule
                  sf_peaks(IncSchemeID,pathwaybuff{1,iSpace},Exp.dt,1,1,Ea,Eb,G1,D1);
                  sf_peaks(IncSchemeID,pathwaybuff{2,iSpace},Exp.dt,2,2,Ea,Eb,G2,D2);
                else
                  sf_peaks(IncSchemeID,buff,Exp.dt,1,1,Ea,Eb,G1,D1);
                  sf_peaks(IncSchemeID,buff,Exp.dt,2,2,Ea,Eb,G2,D2);
                end
              end

            % 2pESEEM ---------------------------------------------------------
            case 1

              if ~all(idealPulse)
                error('Pre-defined 2pESEEM not supported for real pulses.');
              end
              G = prefactor*M;
              D = M;
              T1left = Mt;
              T1right = Mt;
              if (Exp.tau>0)
                Qab = exp(-2i*pi*Exp.tau*Ea)*exp(-2i*pi*Exp.tau*Eb)';
                G = Qab.*G;
                D = conj(Qab).*D;
              end
              
              if Opt.TimeDomain
                if Opt.ProductRule
                  pathwaytd{1,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[1 2],[2 1],Ea,Eb,G,D,T1left,T1right);
                else
                  totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[1 2],[2 1],Ea,Eb,G,D,T1left,T1right);
                end
              else
                if Opt.ProductRule
                  sf_peaks(IncSchemeID,pathwaybuff{1,iSpace},Exp.dt,[1 2],[2 1],Ea,Eb,G,D,T1left,T1right);
                else
                  sf_peaks(IncSchemeID,buff,Exp.dt,[1 2],[2 1],Ea,Eb,G,D,T1left,T1right);
                end
              end

          end % ExperimentID switchyard

        end % if PredefinedExperiment
      end % offset loop
    end % nuclear subspace loop

       
    % Apply product rule
    %--------------------------------------------------------------
    if ~isENDOR
      if Opt.ProductRule

        for iPathway = 1:nPathways
          td_ = 1;
          for iSpace = 1:nSubSpaces
            if Opt.TimeDomain
              thistd = pathwaytd{iPathway,iSpace};
            else
              if (nDimensions==1)
                thistd = ifft(pathwaybuff{iPathway,iSpace})*nPointsF;
                thistd = thistd(1:Exp.nPoints);
              else
                thistd = ifft2dpartial(pathwaybuff{iPathway,iSpace},Exp.nPoints,Opt.PartialIFFT);
                %thistd = thistd;
              end
            end
            td_ = td_.*thistd;
          end
          totaltd = totaltd + Weights(iOri)*OriSelWeight(iT)*td_;
        end

      else

        % TD: nothing
        % FD: IFFT is done after orientation loop

      end % if Opt.ProductRule
    else
      % ENDOR: no action
    end % if ~isENDOR
    
  end % electronic transition loop

end % orientation loop

logmsg(1,'end of orientation/transition loop');
logmsg(1,'%d of %d orientations skipped',nSkippedOrientations,nOrientations);
%=================================================================



%=================================================================
% Postprocessing
%=================================================================
if isENDOR
  
  %endorspc = real(endoroffset - endorspc);
  endorspc = real(endorspc);

  % Normalize modulation signal
  % No need to normalize out the experiment prefactors (due to
  % pulse transfer amplitudes), since they were not included above.
  EqDensityTrace = prod(2*Sys.I+1);
  endorspc = endorspc/(nPathways*EqDensityTrace);
  
  %endorspc = convspec(endorspc,rf(2)-rf(1),Sys.lwEndor);

else
  
  if Opt.ProductRule
    td = totaltd;
  else
    if Opt.TimeDomain
      td = totaltd;
    else
      logmsg(1,'Postprocessing...');
      if (nDimensions==1)
        td = ifft(buff)*numel(buff);
        td = td(1:Exp.nPoints);
      else
        td = ifft2dpartial(buff,Exp.nPoints,Opt.PartialIFFT)*numel(buff);
      end
    end
  end
  
  % Normalize modulation signal
  % No need to normalize out the experiment prefactors (due to
  % pulse transfer amplitudes), since they were not included above.
  EqDensityTrace = prod(2*Sys.I+1);
  td = td/(nPathways*EqDensityTrace);
  
  if ~PredefinedExperiment
    if nPathways>0
      td = td/max(abs(pulseprefactor));
    end
  end

end


%=================================================================
% Time axes, relaxation
%=================================================================
if ~isENDOR
  switch ExperimentID
    case -1
      switch (nDimensions)
        case 1
          t1 = (0:Exp.nPoints-1)*Exp.dt;
        case 2,
          t1 = (0:Exp.nPoints(1)-1)*Exp.dt(1);
          t2 = (0:Exp.nPoints(2)-1)*Exp.dt(2);          
      end
    case 1 % 2p-ESEEM
      t1 = (0:Exp.nPoints-1)*Exp.dt + Exp.tau(1);
    case 2 % 3p-ESEEM
      t1 = (0:Exp.nPoints-1)*Exp.dt + Exp.tau(1) + Exp.T;
    case 3 % 4p-ESEEM
      t1 = (0:Exp.nPoints-1)*Exp.dt + Exp.T;
    case 4 % HYSCORE
      t1 = (0:Exp.nPoints(1)-1)*Exp.dt(1) + Exp.t1;
      t2 = (0:Exp.nPoints(2)-1)*Exp.dt(2) + Exp.t2;
      t = {t1,t2};
    case 5 % Mims ENDOR
      t1 = (0:Exp.nPoints-1)*Exp.dt;
  end

  DecayAdded = 0;
  if any(~isinf(Exp.T1T2))
    logmsg(1,'Adding relaxation decays...');
    T1 = Exp.T1T2(1);
    T2 = Exp.T1T2(2);
    tdecay = [];
    switch ExperimentID
      case 1 % two-pulse ESEEM
        if ~isinf(T2)
          tdecay = exp(-2*t1/T2);
          td = td.*tdecay;
        end
      case 2 % three-pulse ESEEM
        if ~isinf(T1)
          tdecay = exp(-2*Exp.tau/T2)*exp(-t1/T1);
          td = td.*tdecay;
        end
      case 3 % HYSCORE
        if ~isinf(T1)
          tdecay = exp(-t1/T1)'*exp(-t2/T1);
          tdecay = exp(-2*Exp.tau/T2)*tdecay;
          td = td.*tdecay;
        end
    end
    DecayAdded = ~isempty(tdecay);
  end
  
end


%===============================================================
% TD data processing
%===============================================================
logmsg(1,'-final-------------------------------------------------');
logmsg(1,'Data processing...');
if ~isENDOR
  switch nDimensions
    case 1
      if DataProcessing

        % Decay correction
        if (DecayAdded)
          [kk,cc,tdfit] = exponfit(t1,td,1);
          tdx = td - tdfit;
        else
          tdx = td;
        end
        
        % Baseline correction
        tdx = tdx - mean(tdx);
        %plot(t1,tdx); pause;
        %tdx(1) = tdx(1)/2;

        % Apodization
        win = apowin(Opt.Window,numel(tdx)).';
        tdx = tdx.*win;

        % Fourier transformation
        fd = fft(tdx,Opt.ZeroFillFactor*numel(tdx));
        fd = fftshift(fd);
        f1 = fdaxis(Exp.dt,length(fd));
      end


    case 2
      if DataProcessing
        if (DecayAdded)
          tdx = basecorr(td,[1 2],[2 2]);
        else
          tdx = basecorr(td,[1 2],[0 0]);
        end

        w1 = apowin(Opt.Window,Exp.nPoints(1));
        w2 = apowin(Opt.Window,Exp.nPoints(2));

        fd = abs(fftshift(fftn(tdx.*(w1*w2.'),Opt.ZeroFillFactor*Exp.nPoints)));
        f1 = fdaxis(Exp.dt(1),size(fd,1));
        f2 = fdaxis(Exp.dt(2),size(fd,2));
      end

  end

  if max(abs(fd))<1e-300;
    fd = fd*0;
  end

  % Collect output structure
  if DataProcessing
    if (nDimensions==2), out.f1 = f1; out.f2 = f2; else out.f = f1; end
    out.fd = fd;
  else
    out = [];
  end
  
else
  % Collect output structure
  f1 = rf;
  fd = endorspc;
  if max(abs(fd))<1e-300;
    fd = fd*0;
  end
  out = [];
end


%===============================================================
% Output
%===============================================================
if isENDOR
  switch nargout
    case 1, varargout = {fd};
    case 2, varargout = {f1,fd};
    case 3, varargout = {f1,fd,out};
  end
else
  switch nargout
    case 1, varargout = {td};
    case 2, varargout = {t1,td};
    case 3, if (nDimensions==2), varargout = {t1,t2,td}; else varargout = {t1,td,out}; end
    case 4, if (nDimensions==2), varargout = {t1,t2,td,out}; end
  end
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

%==========================================================================
%==========================================================================
%==========================================================================
function Signal = sf_evolve(IncSchemeID,nPoints,dt,incL,incR,Ea,Eb,G,D,T1l,T1r,T2l,T2r,T3l,T3r)

if nargin-9~=2*(numel(incL)-1)
  error('Time-domain evolution kernel: Wrong number of transfer matrices.');
end

NN = numel(G);
Density = G;
Detector = reshape(D.',1,NN);
E = {Ea,Eb};

% Pre-allocate signal array
if numel(nPoints)==1
  Signal = zeros(1,nPoints);
else
  Signal = zeros(nPoints);
end

switch IncSchemeID
  
  case 1 % IncScheme [1], three-pulse ESEEM, five-pulse ESEEM
    FinalDensity = Density(:);
    UUt = exp(-2i*pi*E{incL(1)}*dt)*exp(-2i*pi*E{incR(1)}*dt)';
    UUt = UUt(:);
    for k1 = 1:nPoints
      Signal(k1) = Detector*FinalDensity;
      FinalDensity = UUt.*FinalDensity;
    end

  case 2 % IncScheme [1 1], two-pulse ESEEM, 1D-CP, refocused two-pulse ESEEM
    UUleft  = exp(-2i*pi*E{incL(2)}*dt)*exp(-2i*pi*E{incL(1)}*dt).';
    UUright = exp(+2i*pi*E{incR(1)}*dt)*exp(+2i*pi*E{incR(2)}*dt).';
    for k = 1:nPoints
      % compute density right before detection
      FinalDensity = T1l*Density*T1r;
      % compute trace(Detector*FinalDensity)
      Signal(k) = Detector*FinalDensity(:);
      % Evolve transfer matrices forward and backwards at the same time.
      T1l = UUleft.*T1l; % equivalent to U*T1Left*U'
      T1r = UUright.*T1r; % equivalent to U*T1Right*U'
    end
  
  case 3 % IncScheme [1 -1]
    error('Time-domain evolution for incrementation scheme [1 -1] not supported.');
    
  case 11 % IncScheme [1 2], HYSCORE
    UUt1 = exp(-2i*pi*E{incL(1)}*dt(1))*exp(-2i*pi*E{incR(1)}*dt(1))';
    UUt2 = exp(-2i*pi*E{incL(2)}*dt(2))*exp(-2i*pi*E{incR(2)}*dt(2))';
    UUt2 = UUt2(:);
    for k1 = 1:nPoints(1)
      FinalDensity = reshape(T1l*Density*T1r,NN,1);
      for k2 = 1:nPoints(2)
        Signal(k1,k2) = Detector*FinalDensity;
        FinalDensity = UUt2.*FinalDensity;
      end
      Density = UUt1.*Density;
    end
    
  case 12 % IncScheme [1 2 2]
    error('Time-domain evolution for this experiment not supported.');

  case 13 % IncScheme [1 2 1] 2D-3p
    error('Time-domain evolution for this experiment not supported.');
  
  case 14 % IncScheme [1 1 2], CF-NF
    error('Time-domain evolution for this experiment not supported.');

  case 15 % IncScheme [1 2 2 1], 2D-CP
    error('Time-domain evolution for this experiment not supported.');

  case 16 % IncScheme [1 2 -2 1]
    error('Time-domain evolution for this experiment not supported.');

  case 17 % IncScheme [1 1 2 2], 2D refocused primary ESEEM
    error('Time-domain evolution for this experiment not supported.');

  otherwise
    
    error('Time-domain evolution for this experiment not supported.');
    
end

return





function [FreeL,FreeR,PulsL,PulsR] = pathwayparser(ctp)

%----------------------------------------------
% Coherence transfer pathway parser
%----------------------------------------------

% ctp: a 2D char array with one CTP per row
%   two-pulse ESEEM: '+-'
%   three-pulse ESEEM: ['+a-';'+b-']
%   HYSCORE: ['+ab-';'+ba'] or ['+ab-';'+ba-';'+aa-';'+bb-']
%   etc.

% General superpropagator
%      Sa  Sb  S+  S-
% Sa   aa  ++  a+  +a
% Sb   --  bb  -b  b-
% S+   a-  +b  ab  +-
% S-   -a  b+  -+  ba
%
% E.g. -b in 2nd row and 3rd column describes transfer from
% S+ to Sb. -b is an abbreviation for U- * ... * Ub'.

Display = 0;
idealPulse = (nargout<=2);

% Pre- and post-multiplications
PropL = [1 3 1 3; 4 2 4 2; 1 3 1 3; 4 2 4 2];
PropR = [1 3 3 1; 4 2 2 4; 4 2 2 4; 1 3 3 1];

% Convert ctp string with ab+- to number array with 1234
% General code:
% alpha = a = 1, beta  = b = 2, plus  = + = 3, minus = - = 4
if ischar(ctp)
  code0('ab+-') = [1 2 3 4];
  ctp = code0(ctp);
end

% Get indices for pre- and post-multiplication for pulses
if ~idealPulse
  before = [ones(size(ctp,1),1)*2 ctp(:,1:end-1)];
  idx = ctp + (before-1)*4;
  PulsL = PropL(idx);
  PulsR = PropR(idx);
end

% Get indices for pre- and post-multiplications for free evolutions
idx = 5*ctp - 4;
FreeL = PropL(idx);
FreeR = PropR(idx);

% Display formula
if (Display)
  decode = 'ab+-';
  ecode = 'AB';
  for i = 1:size(ctp,1)
    formula = '';
    for k = 1:size(ctp,2)
      if ~Opt.idealPulse
        formula = ['P' decode(PulsL(i,k)) ' ' formula ' P' decode(PulsR(i,k)) ''''];
      end
      formula = [ecode(FreeL(i,k)) ' ' formula ' ' ecode(FreeR(i,k)) ''''];
    end
    fprintf('%s :   trace(%s)\n',decode(ctp(i,:)),formula);
  end
end
