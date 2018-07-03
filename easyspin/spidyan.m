% spidyan    Simulate spindyanmics during pulse EPR experiments
%
%     [t,Signal] = spidyan(Sys,Exp,Opt)
%     [t,Signal,out] = spidyan(Sys,Exp,Opt)
%
%     Inputs:
%       Sys   ... spin system with electron spin and ESEEM nuclei
%       Exp   ... experimental parameters (time unit us)
%       Opt   ... simulation options
%
%     Outputs:
%       t                 ... time axis
%       Signal            ... simulated signals of detected events
%       out               ... output structure with fields:
%         .FinalState        ... density matrix/matrices at the end of the
%                                experiment
%         .StateTrajectories ... cell array with density matrices from each
%                                timestep during evolution
%         .Events            ... structure containing events

function varargout = spidyan(Sys,Exp,Opt) 

if (nargin==0), help(mfilename); return; end

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
if (nargout>3), error('Too many output arguments.'); end

% Initialize options structure to zero if not given.
if (nargin<3), Opt = struct('unused',NaN); end
if isempty(Opt), Opt = struct('unused',NaN); end

if ~isstruct(Sys)
  error('First input argument (Sys) must be a structure!');
end
if ~isstruct(Exp)
  error('Second input argument (Exp) must be a structure!');
end
if ~isstruct(Opt)
  error('Third input argument (Opt) must be a structure!');
end
% 

% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

logmsg(1,['=begin=spidyan====' datestr(now) '=================']);
logmsg(2,'  log level %d',EasySpinLogLevel);
%----------------------------------------------------------------------
% Preprocess Input
%----------------------------------------------------------------------

% Check for magnetic field
if ~isfield(Exp,'Field')
  logmsg(1,'  no Exp.Field given');
  logmsg(1,'  checking if Sys allows for running without a field...');
  if isfield(Sys,'g')
    error('Exp.Field is required for Sys.g.')
  elseif Sys.ZeemanFreq % && (~isfield(Sys,'S') || all(Sys.S==1/2))
    % if Sys.ZeemanFreq is used, it Exp.Field is not required, however the
    % program needs a reference field to correctly shift frequencies and
    % to calculate nuclear Larmor frequencies (if nuclei are given). The
    % program does so, by "guessing" a field from Sys.ZeemanFreq and gfree
    Exp.Field = (Sys.ZeemanFreq(1)*1e6)*planck/bmagn/(gfree*1e-6);
    logmsg(1,'  Exp.Field is needed for calculating simulation frame frequencies');
    logmsg(1,'  setting Exp.Field to %.1f mT',Exp.Field);
  else
    error('Exp.Field is required for given spin system.')
  end  
elseif length(Exp.Field) ~= 1
  error('Exp.Field must be a single number (in mT).');
end

% Build the Event and Vary structure
Opt.SimulationMode = 'thyme'; % to force s_sequencer to setup for step wise

[Events, Vary, Opt] = s_sequencer(Exp,Opt);

% Adapt Zeeman frequencies for the selected simulation frame
logmsg(1,'-validating spin system--------------------------------');
if isfield(Sys,'ZeemanFreq') && Opt.FrameShift ~= 0
  logmsg(1,'adapting Sys.ZeemanFreq to the simulation frame');
  Sys.ZeemanFreq =  Sys.ZeemanFreq - Opt.FrameShift;
end

% Adapt g-values for the selected simulation frame
if isfield(Sys,'g')
  if isfield(Sys,'S')
    nElectrons = length(Sys.S);
  else
    nElectrons = 1;
  end
  
  if Opt.FrameShift ~= 0
    logmsg(1,'  adapting Sys.g to the simulation frame');
  end
  gshift = (Opt.FrameShift*1e9)*planck/bmagn/(Exp.Field(end)*1e-3);
  
  issize = @(A,siz) all(size(A)==siz);
  fullg = issize(Sys.g,[3*nElectrons 3]);
  if fullg
    gshiftmat = repmat(gshift*eye(3),[nElectrons,1]);
    Sys.g = Sys.g - gshiftmat;
  else
    Sys.g = Sys.g - gshift;
  end
end

% Translate Frequency to g values
if isfield(Sys,'ZeemanFreq')
  logmsg(1,'  translating Sys.ZeemanFreq into g values');
  Sys.ZeemanFreq = Sys.ZeemanFreq*1000; % GHz -> MHz
  if isfield(Sys,'g')
    [~, dgTensor] = size(Sys.g);
  else
    dgTensor = 1; 
  end
  % Recalculates the g value
  for ieSpin = 1 : length(Sys.ZeemanFreq)
    if Sys.ZeemanFreq(ieSpin) ~= 0
      g = (Sys.ZeemanFreq(ieSpin)*1e6)*planck/bmagn/(Exp.Field*1e-3);
      Sys.g(ieSpin,1:dgTensor) = g;
    end
  end
end

% Remove field ZeemanFreq if given, which is spidyan specific
if isfield(Sys,'ZeemanFreq')
  logmsg(1,'  removing Sys.ZeemanFreq from Sys structure');
  Sys = rmfield(Sys,'ZeemanFreq');
end

% Transfer DetOperator from Exp to Opt, to make it available to
% s_propagationsetup
if isfield(Exp,'DetOperator')
  % catch the case when Exp.DetOperator = 'z1' instead of {'z1'}
  if ~iscell(Exp.DetOperator)
    Exp.DetOperator = {Exp.DetOperator};
  end
  Opt.DetOperator = Exp.DetOperator;
end

% Validate and build spin system as well as excitation operators
logmsg(1,'  parsing the Sys structure...');

[Sys, Sigma, DetOps, Events, Relaxation] = s_propagationsetup(Sys,Events,Opt);

% Get Hamiltonian
logmsg(1,'  computing lab frame Hamiltonian');
Ham = sham(Sys,Exp.Field*[0 0 1]);

%----------------------------------------------------------------------
% Propagation
%----------------------------------------------------------------------

% Calls the actual propagation engine
logmsg(1,'-starting propagation----------------------------------');
logmsg(1,'  depending on the Exp setup this may take a while...');
[TimeAxis, RawSignal, FinalState, StateTrajectories, Events] = ...
  s_thyme(Sigma, Ham, DetOps, Events, Relaxation, Vary);
logmsg(1,'  propagation finished!');

%----------------------------------------------------------------------
% Signal Processing
%----------------------------------------------------------------------
logmsg(1,'-validating and processing outout----------------------');
% Signal postprocessing, such as down conversion and filtering and
% checking output of the timeaxis 

nDetOps = numel(DetOps); % number of detection operators

if Opt.SinglePointDetection
  logmsg(1,'  single point detection...');
  if isfield(Exp,'nPoints') && nDetOps ~= 1
    DimSignal =  ndims(RawSignal);
    Signal = permute(RawSignal,[1:(DimSignal-1) DimSignal+1 DimSignal]);
  else
    Signal = RawSignal;
  end
else
  if ~isempty(RawSignal)
    % Check if Exp.DetFreq and Exp.DetOperator have the same length
    if isfield(Exp,'DetFreq')
      logmsg(1,'  veryfing Exp.DetFreq...');
      if isfield(Exp,'DetOperator')
        if length(Exp.DetFreq) ~= length(Exp.DetOperator)
          error('Exp.DetOperator and Exp.DetFreq must contain the same number of elements.')
        end
      elseif length(Exp.DetFreq) ~= 1
        error('If Exp.DetOperator is not provided the default detection operator (S+ for all electrons) is assumed. In this case Exp.DetFreq must contain only one frequency.')
      end
    else
      if ~isfield(Exp,'DetOperator') && isfield(Exp,'mwFreq')
        logmsg(1,'  using Exp.mwFreq as detection frequency');
        Exp.DetFreq = Exp.mwFreq;
      end
    end
    
    logmsg(1,'  processing transients...');
    
    % Adapt FreqTranslation if needed
    FreqTranslation = zeros(1,nDetOps);
    
    if isfield(Exp,'DetFreq') && ~isempty(Exp.DetFreq)
      
      % This adapts the values for DetFreq to simulation frame
      Exp.DetFreq(Exp.DetFreq > 0) = Exp.DetFreq(Exp.DetFreq > 0) - Events{1}.FrameShift;
      Exp.DetFreq(Exp.DetFreq < 0) = Exp.DetFreq(Exp.DetFreq < 0) + Events{1}.FrameShift;
      
      % And then writes them
      FreqTranslation(1:length(Exp.DetFreq)) = - Exp.DetFreq; % To make it a down conversion for negative frequencies add the neg sign
      
    end
    
    % Downconversion/processing of signal
    Signal = signalprocessing(TimeAxis,RawSignal,FreqTranslation);
    
    % If time axis is the same for each data point, it is reduced to a single
    % vector at this point - helps with plotting
    if ~iscell(TimeAxis)
      SizeT = size(TimeAxis);
      linearTimeAxis = reshape(TimeAxis,[prod(SizeT(1:end-1)) SizeT(end)]);
      if size(unique(linearTimeAxis,'rows'),1) == 1
        TimeAxis = linearTimeAxis(1,:);
      end
    end
  else
    logmsg(1,'  nothing was detected...');
    Signal = [];
  end
end

% Reduce the dimensionality of the final state for simulations with only
% one acquisition point
if ndims(FinalState) == 3 && size(FinalState,1) == 1
  FinalState = squeeze(FinalState);
end

% Same for the StateTrajectories (if any exist). If only one vector of
% StateTrajectories was recorded, the nested structure is removed
if ~isempty(StateTrajectories)
  SizeStateTrajectories = size(StateTrajectories);
  if all(SizeStateTrajectories == 1)
    StateTrajectories = StateTrajectories{1};
  end
end

% Assigning outputs
switch nargout
  case 0
    % no output argument - graphical output
    logmsg(1,'-no output requested------------------------------------');
    logmsg(1,'  switching to graphical output');
    if isempty(Signal)
      disp('Detection was switched off, nothing to display.')
    else
      if ~isfield(Exp,'DetOperator')
        % this field will be used as title in the plotting
        Exp.DetOperator = {'Electron Coherence'};
      end
      
      if isfield(Exp,'nPoints')
        nDataPoints = prod(Exp.nPoints);
      else
        nDataPoints = 1;
      end
      
      % Set up figures
      logmsg(1,'  setting up figures...');
      % only make figures if a) transients were detected or b) single point
      % detection with indirect dimension (but not more than 2) and more
      % than one acquisition point
      if ~Opt.SinglePointDetection || (Opt.SinglePointDetection && isfield(Exp,'nPoints') && length(Exp.nPoints) < 3 && nDataPoints > 1)
        for iDetOp = 1 : nDetOps
          figure(iDetOp)
          clf
          if ischar(Exp.DetOperator{iDetOp})
            % title for detection operators that used the sop syntax
            title(['Signal of ' Exp.DetOperator{iDetOp}])
            ylabel(['<' Exp.DetOperator{iDetOp} '>'])
          else
            % title of detection operators that were provided in matrix form
            title(['Signal detected with detection operator no.' num2str(iDetOp)])
            ylabel('<S>')
          end
          xlabel('t (\mus)')
          hold on
        end
      end
      
      if Opt.SinglePointDetection
        % plotting single point detection
        logmsg(1,'  single point detection:');
        if nDataPoints == 1
          % if only a single datapoint was acquired, there is no point in
          % plotting it, instead the output is displayed in the console
          for iDetOp = 1 : nDetOps
            if ischar(Exp.DetOperator{iDetOp})
              string  = ['  Expectation value of ' Exp.DetOperator{iDetOp} ':   '];
            else
              string = ['  Expectation value of operator no.' num2str(iDetOp) ':   '];
            end
            disp([string num2str(Signal(1,iDetOp))]);
          end
        else
          if length(Exp.nPoints) > 2
            % more than two dimensions can not be displayed in a general
            % way
            logmsg(1,'  more than two indirect dimensions - stopping');
            disp('Unable to display more than two indirect dimensions.')
          elseif length(Exp.nPoints) == 1
            % one dimensional case
            logmsg(1,'  creating plot(s) for one indirect dimension');
            for iDetOp = 1 : nDetOps
              figure(iDetOp)
              if length(Exp.Dim1{1,2}) == 1
                % x-axis and its label in case of linear increments
                xvec = Exp.Dim1{1,2}*(0:Exp.nPoints(1)-1);
                xlabel(['\Delta' Exp.Dim1{1,1}])
              else
                % x-axis and its label in case of user provided increments
                xvec = 1:Exp.nPoints(1);
                xlabel(['Data Points ' Exp.Dim1{1,1}])
              end
              plot(xvec,real(squeeze(Signal(:,1,iDetOp))));
            end
          elseif length(Exp.nPoints) == 2
            % two dimensional case
            logmsg(1,'  creating plot(s) for two indirect dimensions');
            for iDetOp = 1 : nDetOps
              if length(Exp.Dim1{1,2}) == 1
                % y-axis and its label in case of linear increments
                yvec = Exp.Dim1{1,2}*(0:Exp.nPoints(1)-1);
                ylabel(['\Delta' Exp.Dim1{1,1}])
              else
                % y-axis and its label in case of userdefined increments
                yvec = 1:Exp.nPoints(1);
                ylabel(['Data Points ' Exp.Dim1{1,1}])
              end
              if length(Exp.Dim2{1,2}) == 1
                % x-axis and its label in case of linear increments
                xvec = Exp.Dim2{1,2}*(0:Exp.nPoints(2)-1);
                xlabel(['\Delta' Exp.Dim2{1,1}])
              else
                % x-axis and its label in case of userdefined increments
                xvec = 1:Exp.nPoints(2);
                xlabel(['Data Points ' Exp.Dim2{1,1}])
              end
              figure(iDetOp)
              surf(xvec,yvec,real(squeeze(Signal(:,:,1,iDetOp))));
              shading flat
            end
          end
        end
      else
        % plotting transients - this can get a little messy
        logmsg(1,'  plotting transient(s)');
        if nDataPoints == 1
          % plotting a single acquisition point
          for iDetOp = 1 : nDetOps
            figure(iDetOp)
            plot(TimeAxis,real(Signal(:,iDetOp)));
          end
        else
          if ~iscell(Signal)
            % plotting if signals are organized in a numeric array
            SignalSize = size(Signal);
            if length(Exp.DetOperator) == 1
              SignalSize(end+1) = 1;
            end
            % reshape such to three dimensions
            Signal = reshape(Signal,[nDataPoints SignalSize(end-1) SignalSize(end)]);
          end
          % set up time axes
          if ~iscell(TimeAxis) && size(TimeAxis,1) == 1
            IndexVec = ones(1,nDataPoints);
          else
            IndexVec = 1:nDataPoints;
          end
          
          for iDetOp = 1 : nDetOps
            figure(iDetOp)
            for iDataPoint = 1 : nDataPoints
              if iscell(Signal)
                % plotting the cell arrays
                plot(TimeAxis{iDataPoint},real(Signal{iDataPoint}(:,iDetOp)));
              else
                % plotting the reshaped numeric arrays
                plot(TimeAxis(IndexVec(iDataPoint),:),real(squeeze(Signal(iDataPoint,:,iDetOp))));
              end
            end
          end
          
        end
        
      end
    end
  case 1
    varargout = {Signal};
  case 2
    varargout = {TimeAxis,Signal};
  case 3
    out.FinalState = FinalState;
    out.StateTrajectories = StateTrajectories;
    out.Events = Events;
    varargout = {TimeAxis,Signal,out};
  otherwise
    error('Incorrect number of output arguments. 1,2, or 3 expected.');
end

%===============================================================
% Report performance
%===============================================================
[Hours,Minutes,Seconds] = elapsedtime(StartTime,clock);
if (Hours>0)
  msg = sprintf('spidyan took %dh%dm%0.3fs',Hours,Minutes,Seconds);
elseif (Minutes>0)
  msg = sprintf('spidyan took %dm%0.3fs',Minutes,Seconds);
else
  msg = sprintf('spidyan took %0.3fs',Seconds);
end
logmsg(1,msg);

logmsg(1,'=end=spidyan======%s=================\n',datestr(now));
