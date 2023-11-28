% esfit   Least-squares fitting for EPR and other data
%
%   esfit(data,fcn,p0,vary)
%   esfit(data,fcn,p0,lb,ub)
%   esfit(___,FitOpt)
%
%   pfit = esfit(___)
%
% Input:
%     data        experimental data, a vector of data points or a cell
%                   array of datasets for global fitting 
%     fcn         simulation/model function handle (@pepper, @garlic, ...
%                   @salt, @chili, or handle to user-defined function)
%                   a user-defined fcn should take a parameter vector p
%                   and return simulated data datasim: datasim = fcn(p)
%     p0          starting values for parameters
%                   EasySpin-style functions: {Sys0,Exp0} or {Sys0,Exp0,Opt0}
%                   other functions: n-element vector
%     vary        allowed variation of parameters
%                   EasySpin-style functions: {vSys} or {vSys,vExp} or {vSys,vExp,vOpt}
%                   other functions: n-element vector
%     lb          lower bounds of parameters
%                   EasySpin-style functions: {lbSys,lbExp} or {lbSys,lbExp,lbOpt}
%                   other functions: n-element vector
%     ub          upper bounds of parameters
%                   EasySpin-style functions: {ubSys,ubExp} or {ubSys,ubExp,ubOpt}
%                   other functions: n-element vector
%     FitOpt      options for esfit
%        .Method  string containing keywords for
%           -algorithm: 'simplex','levmar','montecarlo','genetic','grid','swarm'
%           -target function: 'fcn', 'int', 'dint', 'diff', 'fft'
%        .AutoScale 'lsq', 'maxabs', 'none'; default 'lsq' for
%                 EasySpin simulation functions, otherwise 'none'
%        .BaseLine 0, 1, 2, 3 or [] (or vector for global fitting with different 
%                 baseline order for different datasets)
%        .OutArg  two numbers [nOut iOut], where nOut is the number of
%                 outputs of the simulation function and iOut is the index
%                 of the output argument to use for fitting
%        .Mask    array of 1 and 0 the same size as data vector
%                 values with mask 0 are excluded from the fit 
%                 (cell array for data input consisting of multiple datasets)
%        .weight  array of weights to use when combining residual vectors
%                 of all datasets for global fitting
% Output:
%     fit           structure with fitting results
%       .pfit       fitted parameter vector (contains only active fitting parameters)
%       .pnames     variable names of the fitted parameters
%       .pfit_full  parameter vector including inactive fitting parameters (in GUI)
%       .argsfit    fitted input arguments (if EasySpin-style)
%       .pstd       standard deviation for all parameters
%       .ci95       95% confidence intervals for all parameters
%       .cov        covariance matrix for all parameters
%       .corr       correlation matrix for all parameters
%       .p_start    starting parameter vector for fit
%       .fitraw     fit, as returned by the simulation/model function
%       .fit        fit, including the fitted scale factor
%       .scale      fitted scale factor
%       .baseline   fitted baseline
%       .mask       mask used for fitting
%       .residuals  residuals
%       .ssr        sum of squared residuals
%       .rmsd       root-mean square deviation between input data and fit
%       .bestfithistory  structure containing a list of fitting parameters
%                        corresponding to progressively improved rmsd
%                        values during fitting process and corresponding
%                        rmsd values, for EasySpin functions, a conversion
%                        function returning the EasySpin input structures
%                        given a selected set of fitting parameters is also
%                        included
%

function result = esfit(data,fcn,p0,varargin)

if nargin==0, help(mfilename); return; end

% Start legacy version if GUI is requested with Matlab <R2021b
if verLessThan('Matlab','9.11') && nargout==0
  warning('Your Matlab version (<R2021b) does not support the newest version of the esfit GUI. Switching to legacy version.')
  esfit_legacy(data,fcn,p0,varargin{:})
  return;
end
% Check for display server used on Linux and display warning
if isunix
  [~,output] = system('echo "$XDG_SESSION_TYPE"');
  if ~strcmp(output(1),'x')
    warning(['There are known issues with new Matlab GUIs on Linux systems with display servers other than Xorg.' ...
             'To avoid issues change the display server to Xorg or switch to the legacy version of the esfit GUI using esfit_legacy().'])
  end
end

% Check expiry date
error(eschecker);

% Parse argument list
switch nargin
  case 4
    varyProvided = true;
    pvary = varargin{1};
    Opt = struct;
  case 5
    if isstruct(varargin{2})
      varyProvided = true;
      pvary = varargin{1};
      Opt = varargin{2};
    else
      varyProvided = false;
      lb = varargin{1};
      ub = varargin{2};
      Opt = struct;
    end
  case 6
    varyProvided = false;
    lb = varargin{1};
    ub = varargin{2};
    Opt = varargin{3};
  otherwise
    str = ['You provided %d inputs, but esfit requires 4, 5, or 6.\n',...
           'Examples:\n',...
           '   esfit(data, fcn, p0, pvary)\n',...
           '   esfit(data, fcn, p0, pvary, Opt)\n',...
           '   esfit(data, fcn, p0, lb, ub)\n',...
           '   esfit(data, fcn, p0, lb, ub, Opt)'
           ];
    error(str,nargin);
end

if isempty(Opt)
  Opt = struct;
else
  if ~isstruct(Opt)
    error('Opt (last input argument) must be a structure.');
  end
end

% Set up global structure for data sharing among local functions
global esfitdata    %#ok<*GVMIS>
esfitdata = struct;  % initialize; removes esfitdata from previous esfit run
esfitdata.currFitSet = [];
esfitdata.UserCommand = 0;

% Load utility functions
argspar = esfit_argsparams();


% Experimental data
%-------------------------------------------------------------------------------
if ~iscell(data)
  data = {data};
end
data_vec = [];
for i = 1:numel(data)
  if ~isnumeric(data{i}) || ~isvector(data{i}) || isempty(data{i})
    error('First input must be numeric experimental data in the form of a vector or a cell array of vectors.');
  end
  data_vec = [data_vec; data{i}(:)];
  datasize(i) = numel(data{i});
end

idx = cell(1,numel(data));
idx_ = [0 cumsum(datasize)];
for i = 1:numel(data)
  idx{i} = false(size(data_vec));
  idx{i}(idx_(i)+1:idx_(i+1)) = true;
end

esfitdata.data = data_vec;
esfitdata.idx = idx;
esfitdata.nDataSets = numel(data);
esfitdata.datasize = datasize;

% Model function
%-------------------------------------------------------------------------------
if ~isa(fcn,'function_handle')
  str = 'The simulation/model function (2nd input) must be a function handle.';
  if ischar(fcn)
    error('%s\nUse esfit(data,@%s,...) instead of esfit(data,''%s'',...).',str,fcn,fcn);
  else
    error('%s\nFor example, to use the function pepper(...), use esfit(data,@pepper,...).',str);
  end
end
try
  nargin(fcn);
catch
  error('The simulation/model function given as second input cannot be found. Check the name.');
end

esfitdata.fcn = fcn;
esfitdata.fcnName = func2str(fcn);

esfitdata.lastSetID = 0;

% Determine if the model function is an EasySpin simulation function that
% takes structure inputs
EasySpinFunction = any(strcmp(esfitdata.fcnName,{'pepper','garlic','chili','salt'}));


% Parameters
%-------------------------------------------------------------------------------
structureInputs = isstruct(p0) || iscell(p0);
esfitdata.structureInputs = structureInputs;

% Determine parameter intervals, either from p0 and pvary, or from lower/upper bounds
if structureInputs
  argspar.validargs(p0);
  esfitdata.nSystems = numel(p0{1});
  if varyProvided
    % use p0 and pvary to determine lower and upper bounds
    pinfo = argspar.getparaminfo(pvary);
    argspar.checkparcompatibility(pinfo,p0);
    pvec_0 = argspar.getparamvalues(p0,pinfo);
    pvec_vary = argspar.getparamvalues(pvary,pinfo);
    pvec_lb = pvec_0 - pvec_vary;
    pvec_ub = pvec_0 + pvec_vary;
    % Set lower bound to zero for parameters not allowing negative values
    nonnegFieldNames = {'lw','lwpp','weight','initState'};
    ind = ismember({pinfo.FieldName}.',nonnegFieldNames) & pvec_lb<0;
    pvec_lb(ind) = 0;    
  else
    % use provided lower and upper bounds
    pinfo = argspar.getparaminfo(lb);
    argspar.checkparcompatibility(pinfo,p0);
    argspar.checkparcompatibility(pinfo,ub);
    pvec_0 = argspar.getparamvalues(p0,pinfo);
    pvec_lb = argspar.getparamvalues(lb,pinfo);
    pvec_ub = argspar.getparamvalues(ub,pinfo);
  end
else
  if varyProvided
    pvec_0 = p0;
    pvec_vary = pvary;
    pvec_lb = pvec_0 - pvec_vary;
    pvec_ub = pvec_0 + pvec_vary;
  else
    pvec_0 = p0;
    pvec_lb = lb;
    pvec_ub = ub;
  end
  % Generate parameter names
  for k = numel(p0):-1:1
    pinfo(k).Name = sprintf('p(%d)',k);
  end
end

% Convert all parameter vectors to column vectors
pvec_0 = pvec_0(:);
pvec_lb = pvec_lb(:);
pvec_ub = pvec_ub(:);

% Assert parameter vectors and parameters bounds are valid
nParams = numel(pvec_0);
if numel(pvec_lb)~=nParams
  error('Vector of lower bounds has %d elements, but %d are expected.',numel(pvec_lb),nParams);
end
if numel(pvec_ub)~=nParams
  error('Vector of upper bounds has %d elements, but %d are expected.',numel(pvec_lb),nParams);
end
idx = pvec_lb>pvec_ub;
if any(idx)
  error('Parameter #%d: upper bound cannot be smaller than lower bound.',find(idx,1));
end
idx = pvec_0<pvec_lb;
if any(idx)
  error('Parameter #%d: start value is smaller than lower bound.',find(idx,1));
end
idx = pvec_0>pvec_ub;
if any(idx)
  error('Parameter #%d: start value is larger than upper bound.',find(idx,1));
end

% Eliminate fixed parameters
keep = pvec_lb~=pvec_ub;
pinfo = pinfo(keep);
pvec_0 = pvec_0(keep);
pvec_lb = pvec_lb(keep);
pvec_ub = pvec_ub(keep);
nParameters = numel(pvec_0);

if nParameters==0
  error('No variable parameters to fit.');
end

% Store parameter information in global data structure
esfitdata.args = p0;
esfitdata.pinfo = pinfo;
esfitdata.pvec_0 = pvec_0;
esfitdata.p_start = pvec_0;
esfitdata.pvec_lb = pvec_lb;
esfitdata.pvec_ub = pvec_ub;
esfitdata.nParameters = nParameters;

% Initialize parameter fixing mask used in GUI table
esfitdata.fixedParams = false(1,numel(pvec_0));

% Experimental parameters (for EasySpin functions)
%-------------------------------------------------------------------------------
if EasySpinFunction

  if numel(data)>1
    error('Cannot use EasySpin functions for global fitting directly. Write a custom model function.');
  end

  if ~iscell(p0) || numel(p0)<2
    error('The third input must contain the initial parameters, e.g. {Sys0,Exp} or {Sys0,Exp,Opt}.');
  end

  % Check or set Exp.nPoints
  if isfield(p0{2},'nPoints')
    if p0{2}.nPoints~=numel(data{1})
      error('Exp.nPoints is %d, but the data vector has %d elements.',...
        p0{2}.nPoints,numel(data{1}));
    end
  else
    p0{2}.nPoints = numel(data{1});
  end
  
  % For field and frequency sweeps, require manual field range (to prevent
  % users from comparing sim and exp spectra with different ranges)
  if ~any(isfield(p0{2},{'Range','CenterSweep','mwRange','mwCenterSweep'}))
    error('Please specify field or frequency range, in Exp.Range/Exp.mwRange or in Exp.CenterSweep/Exp.mwCenterSweep.');
  end

  % Get x-axis for plotting in GUI
  if isfield(p0{2},'mwRange') || isfield(p0{2},'mwCenterSweep')
    rangefield = 'mwRange';
    centersweepfield = 'mwCenterSweep';
  elseif isfield(p0{2},'Range') || isfield(p0{2},'CenterSweep')
    rangefield = 'Range';
    centersweepfield = 'CenterSweep';
  end
  if isfield(p0{2},centersweepfield) && all(~isnan(p0{2}.(centersweepfield)))
    range = p0{2}.(centersweepfield)(1) + [-1 1]/2*p0{2}.(centersweepfield)(2);
  else
    range = p0{2}.(rangefield);
  end
  if any(range<0) || diff(range)<=0 || any(~isfinite(range)) || any(~isreal(range))
    error('Invalid sweep range! Check Exp.(mw)CenterSweep or Exp.(mw)Range.');
  end
  Opt.x = linspace(range(1),range(2),p0{2}.nPoints);

end

if structureInputs
  esfitdata.p2args = @(pars) argspar.setparamvalues(p0,pinfo,pars);
end


% Options
%===============================================================================
if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 1;
end

if isfield(Opt,'Scaling')
  error('Fitting option Opt.Scaling has been replaced by Opt.AutoScale.');
end

if ~isfield(Opt,'OutArg')
  esfitdata.nOutArguments = abs(nargout(esfitdata.fcn));
  esfitdata.OutArgument = esfitdata.nOutArguments;
else
  if numel(Opt.OutArg)~=2
    error('Opt.OutArg must contain two values [nOut iOut].');
  end
  if Opt.OutArg(2)>Opt.OutArg(1)
    error('Opt.OutArg: second number cannot be larger than first one.');
  end
  esfitdata.nOutArguments = Opt.OutArg(1);
  esfitdata.OutArgument = Opt.OutArg(2);  
end

if ~isfield(Opt,'Method')
  Opt.Method = 'simplex fcn';
end
if EasySpinFunction
  if isfield(p0{2},'Harmonic') && p0{2}.Harmonic>0
    Opt.TargetID = 2; % integral
  else
    if strcmp(esfitdata.fcnName,'pepper') || strcmp(esfitdata.fcnName,'garlic')
      Opt.TargetID = 2; % integral
    end
  end
end

keywords = split(Opt.Method,' ');
for k = 1:numel(keywords)
  switch keywords{k}
    case 'simplex',    Opt.AlgorithmID = 1;
    case 'levmar',     Opt.AlgorithmID = 2;
    case 'montecarlo', Opt.AlgorithmID = 3;
    case 'genetic',    Opt.AlgorithmID = 4;
    case 'grid',       Opt.AlgorithmID = 5;
    case 'swarm',      Opt.AlgorithmID = 6;
    case 'lsqnonlin',  Opt.AlgorithmID = 7;
      
    case 'fcn',        Opt.TargetID = 1;
    case 'int',        Opt.TargetID = 2;
    case 'dint',       Opt.TargetID = 3;
    case 'diff',       Opt.TargetID = 4;
    case 'fft',        Opt.TargetID = 5;
    otherwise
      error('Unknown ''%s'' in Opt.Method.',keywords{k});
  end
end

if ~isfield(Opt,'TargetID')
  Opt.TargetID = 1;
end

AlgorithmNames{1} = 'Nelder-Mead simplex';
AlgorithmNames{2} = 'Levenberg-Marquardt';
AlgorithmNames{3} = 'Monte Carlo';
AlgorithmNames{4} = 'genetic algorithm';
AlgorithmNames{5} = 'grid search';
AlgorithmNames{6} = 'particle swarm';
AlgorithmNames{7} = 'lsqnonlin';
esfitdata.AlgorithmNames = AlgorithmNames;

AlgorithmAbbrev{1} = 'simplex';
AlgorithmAbbrev{2} = 'levmar';
AlgorithmAbbrev{3} = 'montecarlo';
AlgorithmAbbrev{4} = 'genetic';
AlgorithmAbbrev{5} = 'grid';
AlgorithmAbbrev{6} = 'swarm';
AlgorithmAbbrev{7} = 'lsqnonlin';
esfitdata.AlgorithmAbbrev = AlgorithmAbbrev;

TargetNames{1} = 'data as is';
TargetNames{2} = 'integral';
TargetNames{3} = 'double integral';
TargetNames{4} = 'derivative';
TargetNames{5} = 'Fourier transform';
esfitdata.TargetNames = TargetNames;

TargetAbbrev{1} = 'fcn';
TargetAbbrev{2} = 'int';
TargetAbbrev{3} = 'dint';
TargetAbbrev{4} = 'diff';
TargetAbbrev{5} = 'fft';
esfitdata.TargetAbbrev = TargetAbbrev;

% Mask
if ~isfield(Opt,'Mask')
  Opt.Mask = true(size(data_vec));
else
  if iscell(Opt.Mask)
    for i = 1:numel(Opt.Mask)
      Opt.Mask{i} = Opt.Mask{i}(:);
    end
    Opt.Mask = cat(1,Opt.Mask{:});
  end
  Opt.Mask = logical(Opt.Mask(:));
  if numel(Opt.Mask)~=numel(data_vec)
    error('Opt.Mask has %d elements, but the data has %d elements.',numel(Opt.Mask),numel(data));
  end
end
Opt.useMask = true;

% Scale fitting
if ~isfield(Opt,'AutoScale')
  if EasySpinFunction || esfitdata.nDataSets>1
    Opt.AutoScale = 'lsq';
  else
    Opt.AutoScale = 'none';
  end
end
switch Opt.AutoScale
  case 'none', AutoScaleID = 0;
  case 'lsq', AutoScaleID = 1;
  case 'maxabs', AutoScaleID = 2;
  otherwise, error('Unknown setting for Opt.AutoScale - possible values are ''lsq'', ''maxabs'' and ''none''.');
end
Opt.AutoScaleID = AutoScaleID;

esfitdata.AutoScaleSettings = {0, 1, 2};
esfitdata.AutoScaleStrings = {'none', 'lsq', 'maxabs'};

% Baseline correction
if ~isfield(Opt,'BaseLine')
  Opt.BaseLine = [];
end
if isempty(Opt.BaseLine)
  Opt.BaseLine = -ones(1,esfitdata.nDataSets);  
elseif numel(Opt.BaseLine)==1 && esfitdata.nDataSets~=1
  Opt.BaseLine = Opt.BaseLine*ones(1,esfitdata.nDataSets);
else
  Opt.BaseLine = Opt.BaseLine;
end
esfitdata.BaseLine = Opt.BaseLine;
if numel(Opt.BaseLine)~=esfitdata.nDataSets
  error('The number of entries in Opt.BaseLine must be equal to the number of datasets (%d).',esfitdata.nDataSets);
end

esfitdata.BaseLineSettings = {-1, 0, 1, 2, 3};
esfitdata.BaseLineStrings = {'none', 'offset', 'linear', 'quadratic', 'cubic'};

% Uncertainty calculation
if ~isfield(Opt,'CalculateUncertainties')
  Opt.CalculateUncertainties = 1;
end

% Weights for global fitting
if ~isfield(Opt,'weight')
  Opt.weight = ones(1,esfitdata.nDataSets);
end
if numel(Opt.weight)~=esfitdata.nDataSets
  error('The number of elements in Opt.weight must be equal to the number of datasets.');
end

weight = zeros(1,sum(esfitdata.datasize));
for i = 1:numel(data)
  weight(esfitdata.idx{i}) = Opt.weight(i);
end
esfitdata.weight = weight;

% x axis for plotting
showxaxis = false;
if isfield(Opt,'x')
  if iscell(Opt.x)
    x = [];
    for i = 1:numel(Opt.x)
      x = [x; Opt.x{i}(:)];
    end
    Opt.x = x;
  end
  if numel(Opt.x)~=numel(data_vec)
    error('The size of Opt.x must match the size of the experimental data.');
  end
  showxaxis = true;
else
  for i = 1:esfitdata.nDataSets
    Opt.x(esfitdata.idx{i}) = 1:esfitdata.datasize(i);
  end
end
Opt.x = Opt.x(:);

esfitdata.rmsdhistory = [];

esfitdata.besthistory.rmsd = [];
esfitdata.besthistory.par = [];

% Internal parameters
if ~isfield(Opt,'PlotStretchFactor'), Opt.PlotStretchFactor = 0.05; end
if ~isfield(Opt,'maxParameters'), Opt.maxParameters = 30; end

if esfitdata.nParameters>Opt.maxParameters
    error('Cannot fit more than %d parameters simultaneously.',...
        Opt.maxParameters);
end
Opt.IterationPrintFunction = @iterationprint;

% Setup GUI and return if in interactive mode or close GUI
%-------------------------------------------------------------------------------
interactiveMode = nargout==0;
Opt.InfoPrintFunction = @(str) infoprint(str,interactiveMode);
esfitdata.Opts = Opt;
if interactiveMode
  global gui %#ok<TLEV> 
  gui = struct;  % initialize
  gui.showxaxis = showxaxis;
  setupGUI(data);
  return
else
  % Close GUI window if running in script mode
  hFig = findall(0,'Tag','esfitFigure');
  if ~isempty(hFig)
    close(hFig);
    clear hFig
    esfitdata.UserCommand = 0;
  end
  clear global gui
end

% Report parsed inputs
%-------------------------------------------------------------------------------
if esfitdata.Opts.Verbosity>=1
  nDataSets = esfitdata.nDataSets;
  fprintf('-- esfit ------------------------------------------------\n');
  fprintf('Number of datasets:       %d\n',nDataSets);
  for i = 1:nDataSets
    fprintf('Data set %d:           %d points\n',i,esfitdata.datasize(i));
  end
  fprintf('Model function name:      %s\n',esfitdata.fcnName);
  fprintf('Number of fit parameters: %d\n',esfitdata.nParameters);
  fprintf('Minimization algorithm:   %s\n',esfitdata.AlgorithmNames{esfitdata.Opts.AlgorithmID});
  fprintf('Residuals computed from:  %s\n',esfitdata.TargetNames{esfitdata.Opts.TargetID});
  fprintf('Autoscaling:              %s\n',esfitdata.Opts.AutoScale);
  fprintf('---------------------------------------------------------\n');
end

% Run least-squares fitting
%-------------------------------------------------------------------------------
result = runFitting();

clear global esfitdata
clear global gui

end

%===============================================================================
%===============================================================================
%===============================================================================

%===============================================================================
% Run fitting algorithm
%===============================================================================
function result = runFitting(useGUI)

if nargin<1, useGUI = false; end

global esfitdata gui
data_ = esfitdata.data;
fixedParams = esfitdata.fixedParams;
activeParams = ~fixedParams;
Verbosity = esfitdata.Opts.Verbosity;

% Reset best fit history
esfitdata.besthistory.rmsd = [];
esfitdata.besthistory.par = [];

if useGUI
  esfitdata.modelErrorHandler = @(ME) GUIErrorHandler(ME);
else
  esfitdata.modelErrorHandler = @(ME) error('\nThe model simulation function raised the following error:\n  %s\n',ME.message);
end

% Set starting point
%-------------------------------------------------------------------------------
p_start = esfitdata.p_start;
lb = esfitdata.pvec_lb;
ub = esfitdata.pvec_ub;

esfitdata.best.rmsd = inf;
esfitdata.best.rmsdtarget = inf;

% Run minimization over space of active parameters
%-------------------------------------------------------------------------------
fitOpt = esfitdata.Opts;
nActiveParams = sum(activeParams);
if nActiveParams>0
  if Verbosity>=1
    msg = sprintf('Running optimization algorithm with %d active parameters...',nActiveParams);
    if useGUI
      set(gui.nParamsField,'Value',nActiveParams);
      set(gui.nParamsFieldPopup,'Value',nActiveParams);
      updateLogBox(msg)
    else
      disp(msg);
    end
  end
  if useGUI
    fitOpt.IterFcn = @iterupdateGUI;
  end
  fitOpt.track = true;
  if useGUI && (fitOpt.AlgorithmID==6 || fitOpt.AlgorithmID==7)
    iterupdate = true;
  else
    iterupdate = false;
  end
  residualfun = @(x) residuals_(x,fitOpt,iterupdate);
  rmsdfun = @(x) rmsd_(x,fitOpt,iterupdate);
  p0_active = p_start(activeParams);
  lb_active = lb(activeParams);
  ub_active = ub(activeParams);
  switch fitOpt.AlgorithmID
    case 1 % Nelder-Mead simplex
      pfit_active = esfit_simplex(rmsdfun,p0_active,lb_active,ub_active,fitOpt);
    case 2 % Levenberg-Marquardt
      pfit_active = esfit_levmar(residualfun,p0_active,lb_active,ub_active,fitOpt);
    case 3 % Monte Carlo
      pfit_active = esfit_montecarlo(rmsdfun,lb_active,ub_active,fitOpt);
    case 4 % Genetic
      pfit_active = esfit_genetic(rmsdfun,lb_active,ub_active,fitOpt);
      pfit_active = pfit_active(:);
    case 5 % Grid search
      pfit_active = esfit_grid(rmsdfun,lb_active,ub_active,fitOpt);
    case 6 % Particle swarm
      pfit_active = esfit_swarm(rmsdfun,lb_active,ub_active,fitOpt);
    case 7 % lsqnonlin from Optimization Toolbox
      [pfit_active,~,~,~,output] = lsqnonlin(residualfun,p0_active,lb_active,ub_active);
      info.bestx = pfit_active;
      info.newbest = true;
      iterupdateGUI(info);
      if Verbosity>=1 && useGUI && isfield(info,'msg')
        updateLogBox(output.message);
      end
  end
  pfit = p_start;
  pfit(activeParams) = pfit_active;
else
  if Verbosity>=1
    msg = 'No active parameters; skipping optimization';
    if useGUI
      updateLogBox(msg);
    else
      disp(msg);
    end
  end
  pfit = p_start;
end

if isfield(esfitdata,'modelEvalError') && esfitdata.modelEvalError
  return;
end

if esfitdata.structureInputs
  argsfit = esfitdata.p2args(pfit);
else
  argsfit = [];
end

% Get best-fit spectrum
fit = esfitdata.best.fit;  % bestfit is set in residuals_
scale = esfitdata.best.scale;  % bestscale is set in residuals_
for i = 1:esfitdata.nDataSets
  idx = esfitdata.idx{i};
  fitraw(idx,:) = fit(idx)/scale(i);
end
baseline = esfitdata.best.baseline;
baselinetype = esfitdata.best.baselinetype;

% Calculate metrics for goodness of fit
%-------------------------------------------------------------------------------
rmsd0 = esfitdata.best.rmsd;
residuals0 = esfitdata.best.residuals;
ssr0 = sum(abs(residuals0).^2); % sum of squared residuals

% Noise estimated from residuals (assumes excellent fit)
noise_std = std(residuals0);
% Noise estimate following DER_SNR algorithm by Stoehr et al.
N = numel(data_(:));
noise_der = 1.482602/sqrt(6)*median(abs(2.0*data_(3:N-2) - data_(1:N-4) - data_(5:N)));

% Reduced chi square
if noise_std>noise_der
  red_chisquare = (1/N)*sum(residuals0.^2)/noise_der^2;
  chisquareinfo = '(using noise der estimate)';
else
  red_chisquare = (1/N)*sum(residuals0.^2)/noise_std^2;
  chisquareinfo = '(using noise std estimate; upper limit)';
end

% Calculate parameter uncertainties
%-------------------------------------------------------------------------------
calculateUncertainties = esfitdata.Opts.CalculateUncertainties && esfitdata.UserCommand==0 && nActiveParams>0;
if calculateUncertainties
  if Verbosity>=1
    if useGUI
      clear msg
      msg{1} = '';
      msg{2} = 'Calculating parameter uncertainties...';
      msg{3} = '  Estimating Jacobian...';
      updateLogBox(msg);
    else
      disp('Calculating parameter uncertainties...');
      disp('  Estimating Jacobian...');
    end
  end
  %maxRelStep = min((ub-pfit),(pfit-lb))./pfit;
  fitOpt.track = false;
  residualfun = @(x)residuals_(x,fitOpt,useGUI);
  J = jacobianest(residualfun,pfit_active);
  if ~any(isnan(J(:))) && ~isempty(J)
    if Verbosity>=1
      msg = '  Calculating parameter covariance matrix...';
      if useGUI
        updateLogBox(msg);
      else
        disp(msg);
      end
    end

    % Calculate covariance matrix and standard deviations
    residuals = calculateResiduals(fit(:),data_(:),esfitdata.Opts.TargetID,esfitdata.idx);
    residuals = residuals.'; % col -> row
    covmatrix = hccm(J,residuals,'HC1');
    pstd = sqrt(diag(covmatrix));

    % Calculate confidence intervals
    norm_icdf = @(p)-sqrt(2)*erfcinv(2*p); % inverse of standard normal cdf
    ci = @(pctl)norm_icdf(1/2+pctl/2)*sqrt(diag(covmatrix));
    pctl = 0.95;
    ci95 = pfit_active + ci(pctl)*[-1 1];

    % Calculate correlation matrix
    if Verbosity>=1
      msg = '  Calculating parameter correlation matrix...';
      if useGUI
        updateLogBox(msg);
      else
        disp(msg);
      end
    end
    Q = diag(diag(covmatrix).^(-1/2));
    corrmatrix = Q*covmatrix*Q;

    % Report fit results
    %---------------------------------------------------------------------------
    if esfitdata.Opts.Verbosity>=1 || useGUI
      clear msg
      msg{1} = '';
      if esfitdata.Opts.AutoScaleID~=0
        msg{end+1} = ' ';
        msg{end+1} = sprintf('Fitted scale:      %g\n',scale);
      end
      if ~useGUI
        msg{end+1} = 'Parameters:';
        msg{end+1} = printparlist(pfit_active,esfitdata.pinfo,pstd,ci95);
        msg{end+1} = ' ';
      end
      if ~isempty(corrmatrix) && numel(pfit_active)>1
        msg{end+1} = sprintf('Correlation matrix:');
        Sigma = corrmatrix;
        msg{end+1} = sprintf(['    ',repmat('%f  ',1,size(Sigma,1)),'\n'],Sigma);
        triuCorr = triu(abs(Sigma),1);
        msg{end+1} = sprintf('Strongest correlations:');
        [~,idx] = sort(triuCorr(:),'descend');
        [i1,i2] = ind2sub(size(Sigma),idx);
        np = numel(pfit_active);
        parind = 1:numel(activeParams);
        parind = parind(activeParams);
        for k = 1:min(5,(np-1)*np/2)
          msg{end+1} = sprintf('    p(%d)-p(%d):    %g',parind(i1(k)),parind(i2(k)),Sigma(i1(k),i2(k))); %#ok<*AGROW> 
        end
        if any(reshape(triuCorr,1,[])>0.8)
          msg{end+1} = '    WARNING! Strong correlations between parameters.';
        end
      end
      msg{end+1} = ' ';
      msg{end+1} = 'Goodness of fit:';
      msg{end+1} = sprintf('   ssr             %g',ssr0);
      msg{end+1} = sprintf('   rmsd            %g',rmsd0);
      msg{end+1} = sprintf('   noise std       %g (estimated from residuals; assumes excellent fit)',noise_std);
      msg{end+1} = sprintf('   noise der       %g (estimated using der_snr algorithm)',noise_der);
      msg{end+1} = sprintf('   red chi-square  %g %s',red_chisquare,chisquareinfo);
     if useGUI
        msg{end+1} = '';
        updateLogBox(msg);
      else
        disp(repmat('-',1,110));
        for i = 1:numel(msg)
          disp(msg{i});
        end
        disp(repmat('-',1,110));
      end
    end
  else
    if Verbosity>=1
      if isempty(J)
        msg = '  Jacobian estimation interrupted by user, cannot calculate parameter uncertainties.';
      else
        msg = '  NaN elements in Jacobian, cannot calculate parameter uncertainties.';
      end
      if useGUI
        updateLogBox(msg);
      else
        disp(msg);
      end
    end
    pstd = [];
    ci95 = [];
    covmatrix = [];
    corrmatrix = [];
  end
else
  if Verbosity>=1
    if esfitdata.Opts.CalculateUncertainties
      msg = 'Fitting stopped by user. Skipping uncertainty quantification.';
    else
      msg = 'Skipping uncertainty quantification as requested.';
    end    
    if useGUI
      updateLogBox({msg,''});
    else
      disp(msg);
    end
  end
  pstd = [];
  ci95 = [];
  covmatrix = [];
  corrmatrix = [];
end

% Assemble output structure
%-------------------------------------------------------------------------------
result.algorithm = esfitdata.AlgorithmAbbrev{esfitdata.Opts.AlgorithmID};
result.target = esfitdata.TargetAbbrev{esfitdata.Opts.TargetID};
if esfitdata.nDataSets>1
  for k = 1:esfitdata.nDataSets
    idx = esfitdata.idx{k};
    result.fit{k} = fit(idx);
    result.fitraw{k} = fitraw(idx);
    result.mask{k} = esfitdata.Opts.Mask(idx);
    result.residuals{k} = residuals0(idx);
    result.baseline{k} = baseline(idx);
  end
else
  result.fit = fit(:);
  result.fitraw = fitraw(:);
  result.mask = esfitdata.Opts.Mask;
  result.residuals = residuals0(:);
  result.baseline = baseline(:);
end

result.baselinetype = baselinetype;
result.scale = scale;

result.bestfithistory.rmsd = esfitdata.besthistory.rmsd;
result.bestfithistory.pfit = esfitdata.besthistory.par;
if esfitdata.structureInputs
  result.bestfithistory.pfit2structs = esfitdata.p2args;
end

result.pnames = {esfitdata.pinfo.Name}.';
result.pnames = result.pnames(activeParams);
result.p_start = p_start;
result.p_fixed = fixedParams;
result.pfit = pfit_active;
result.pfit_full = pfit;

result.argsfit = argsfit;

result.pstd = pstd;
result.ci95 = ci95;
result.cov = covmatrix;
result.corr = corrmatrix;

result.ssr = ssr0;
result.rmsd = rmsd0;
result.redchisquare = red_chisquare;

esfitdata.best.fit = fit;
fieldnames = {'pstd','ci95','pfit'};
for i = 1:numel(fieldnames)
  esfitdata.best.(fieldnames{i}) = result.(fieldnames{i});
end

end
%===============================================================================

%===============================================================================
function [rmsd,userstop] = rmsd_(x,Opt,iterupdate)
[~,rmsd,userstop] = residuals_(x,Opt,iterupdate);
end
%===============================================================================

%===============================================================================
function [residuals,rmsd,userstop] = residuals_(x,Opt,iterupdate)

global esfitdata

userstop = esfitdata.UserCommand~=0;

expdata = esfitdata.data;

if esfitdata.Opts.useMask
  mask = Opt.Mask;
else
  mask = true(size(Opt.Mask));
end

% Assemble full parameter vector
%-------------------------------------------------------------------------------
par = esfitdata.p_start;
active = ~esfitdata.fixedParams;
par(active) = x;

% Evaluate model function
%-------------------------------------------------------------------------------
out = cell(1,esfitdata.nOutArguments);
try
  if esfitdata.structureInputs
    args = esfitdata.p2args(par);
    [out{:}] = esfitdata.fcn(args{:});
  else
    [out{:}] = esfitdata.fcn(par);
  end
  esfitdata.modelEvalError = false;
catch ME
  esfitdata.modelErrorHandler(ME);
  esfitdata.modelEvalError = true;
  return
end

simdata = out{esfitdata.OutArgument}; % pick appropriate output argument
if ~iscell(simdata)
  simdata = {simdata};
end

if numel(simdata)~=esfitdata.nDataSets
  error('\n  Experimental and model data have unequal number of datasets:\n    experimental: %d\n    model: %d\n',...
    numel(expdata),numel(simdata));
end

simdata_vec = [];
for i = 1:numel(simdata)
  simdata_vec = [simdata_vec; simdata{i}(:)];
end

if numel(simdata_vec)~=numel(expdata)
  error('\n  Experimental data and model have unequal total number of points:\n    experimental: %d\n    model: %d\n',...
    numel(expdata),numel(simdata_vec));
end

% Rescale simulated data if scale should be ignored; include baseline if wanted
%-------------------------------------------------------------------------------
order = Opt.BaseLine;
baseline = zeros(size(simdata_vec));
simscale = ones(1,esfitdata.nDataSets);
for k = 1:esfitdata.nDataSets
  baselinetype(k) = esfitdata.BaseLineStrings([esfitdata.BaseLineSettings{:}]==order(k));
  if order(k)~=-1
    N = esfitdata.datasize(k);
    x = (1:N).'/N;
    D = x.^(0:order(k));  % each column a x^j monomial vector
    idx = esfitdata.idx{k};
    switch Opt.AutoScaleID
      case 0 % 'none'
        coeffs = D(mask(idx),:)\(expdata(mask&idx)-simdata_vec(mask&idx));
        baseline(idx) = D*coeffs;
        simdata_vec(idx) = simdata_vec(idx) + baseline(idx);
      case 1 % 'lsq'
        D = [simdata_vec(idx) D];
        coeffs = D(mask(idx),:)\expdata(mask & idx);
        coeffs(1) = abs(coeffs(1));
        baseline(idx) = D(:,2:end)*coeffs(2:end);
        simdata_vec(idx) = D*coeffs;
        simscale(k) = coeffs(1);
      case 2 % 'maxabs'
        D = [simdata_vec(idx) D];
        coeffs = D(mask(idx),:)\expdata(mask & idx);
        baseline(idx) = D(:,2:end)*coeffs(2:end);
        coeff = max(abs(simdata_vec(mask & idx)))\max(abs(expdata(mask & idx)-baseline(mask & idx)));
        simdata_vec(idx) = coeff*simdata_vec(idx)+baseline(idx);
        simscale(k) = coeff;
    end
  else
    switch Opt.AutoScaleID
      case 1 % 'lsq'
        idx = esfitdata.idx{k};
        coeffs = simdata_vec(mask & idx)\expdata(mask & idx);
        coeffs(1) = abs(coeffs(1));
        simdata_vec(idx) = simdata_vec(idx)*coeffs;
        simscale(k) = coeffs(1);
      case 2 % 'maxabs'
        idx = esfitdata.idx{k};
        coeff = max(abs(simdata_vec(mask & idx)))\max(abs(expdata(mask & idx)));
        simdata_vec(idx) = simdata_vec(idx)*coeff;
        simscale(k) = coeff;    
    end
  end
end

% Compute residuals
%-------------------------------------------------------------------------------
[residuals,residuals0] = calculateResiduals(simdata_vec,expdata,Opt.TargetID,esfitdata.idx,mask(:));
rmsd = sqrt(mean(abs(residuals).^2.*esfitdata.weight(:)));
rmsd0 = sqrt(mean(abs(residuals0).^2));

esfitdata.curr.rmsd = rmsd0;
esfitdata.curr.sim = simdata_vec;
esfitdata.curr.par = par;
esfitdata.curr.scale = simscale;
esfitdata.curr.baseline = baseline;
esfitdata.curr.baselinetype = baselinetype;

% Keep track of errors
%-------------------------------------------------------------------------------
if Opt.track
  esfitdata.rmsdhistory = [esfitdata.rmsdhistory rmsd0];

  isNewBest = rmsd<esfitdata.best.rmsdtarget;
  if isNewBest
    esfitdata.best.residuals = residuals0;
    esfitdata.best.rmsdtarget = rmsd;
    esfitdata.best.rmsd = rmsd0;
    esfitdata.best.fit = simdata_vec;
    esfitdata.best.scale = simscale;
    esfitdata.best.par = par;
    esfitdata.best.baseline = baseline;
    esfitdata.best.baselinetype = baselinetype;
    
    esfitdata.besthistory.rmsd = [esfitdata.besthistory.rmsd rmsd0];
    esfitdata.besthistory.par = [esfitdata.besthistory.par par];
    
  end
  
  if iterupdate
    info.newbest = isNewBest;
    iterupdateGUI(info);
  end

end

end
%===============================================================================

%===============================================================================
function [residuals,residuals0] = calculateResiduals(A,B,mode,idx,includemask)

if nargin>4
  A(~includemask) = 0;
  B(~includemask) = 0;
end

% ignore residual if A or B is NaN
idxNaN = isnan(A) | isnan(B);
A(idxNaN) = 0;
B(idxNaN) = 0;

for i = 1:numel(idx)
  residuals0{i} = A(idx{i}) - B(idx{i});
  switch mode
    case 1  % fcn
      residuals{i} = residuals0{i};
    case 2  % int
      residuals{i} = cumsum(residuals0{i});
    case 3  % iint
      residuals{i} = cumsum(cumsum(residuals0{i}));
    case 4  % fft
      residuals{i} = abs(fft(residuals0{i}));
    case 5  % diff
      residuals{i} = deriv(residuals0{i});
  end
end
residuals0 = cat(1,residuals0{:});
residuals = cat(1,residuals{:});

end
%===============================================================================

%===============================================================================
% Print parameters, and their uncertainties if availabe.
function str = printparlist(par,pinfo,pstd,pci95)

nParams = numel(par);

maxNameLength = max(arrayfun(@(x)length(x.Name),pinfo));
indent = '   ';

printUncertainties = nargin>2 && ~isempty(pstd);
if printUncertainties
  str = [indent sprintf('    name%svalue        standard deviation        95%% confidence interval',repmat(' ',1,max(maxNameLength-4,0)+2))];
  for p = 1:nParams
    pname = pad(pinfo(p).Name,maxNameLength);
    str_ = sprintf('%2.0i  %s  %-#12.7g %-#12.7g (%6.3f %%)   %-#12.7g - %-#12.7g',p,pname,par(p),pstd(p),pstd(p)/par(p)*100,pci95(p,1),pci95(p,2));
    str = [str newline indent str_];
  end
else
  str = [indent sprintf('    name%svalue',repmat(' ',1,max(maxNameLength-4,0)+2))];
  for p = 1:nParams
    pname = pad(pinfo(p).Name,maxNameLength);
    str_ = sprintf('%2.0i  %s  %-#12.7g',pname,par(p));
    str = [str newline indent str_];
  end
end

if nargout==0
  disp(str);
end

end
%===============================================================================

%===============================================================================
function infoprint(str,useGUI)
if useGUI
  updateLogBox(str);
else
  if iscell(str)
    for i = 1:numel(str)
      disp(str{i});
    end
  else
    disp(str);
  end
end
end
%===============================================================================

%===============================================================================
function setupGUI(data)

global esfitdata gui

Opt = esfitdata.Opts;

% Main figure
%-------------------------------------------------------------------------------
hFig = findall(0,'Tag','esfitFigure');
if ~isempty(hFig)
  delete(hFig);
  clear hFig
end
hPopup = findall(0,'Tag','algorithmpopup');
if ~isempty(hPopup)
  delete(hPopup);
  clear hPopup
end
gui.Fig = uifigure('Tag','esfitFigure');
if ~strcmp(gui.Fig.WindowStyle,'normal')
  gui.Fig.WindowStyle = 'normal';
end
set(gui.Fig,'Visible','off')

sz = [1330 800]; % figure size
screensize = get(0,'ScreenSize');
scalefact = min(0.9*(screensize(3:4)/sz));
if scalefact>1
  scalefact = 1;
end
sz = sz*scalefact;
xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the screen vertically
set(gui.Fig,'position',[xpos, ypos, sz(1), sz(2)],'units','pixels');
set(gui.Fig,'MenuBar','none');
% set(gui.Fig,'Resize','off');
set(gui.Fig,'Name','esfit - Least-Squares Fitting','NumberTitle','off');
set(gui.Fig,'CloseRequestFcn',@closeGUI);
  
spacing = 20*scalefact;
hPtop = 190*scalefact;
hElement = 22*scalefact; % height of popup menu, checkboxes, small buttons

fontsizetbl = 11;

% Set up grid layout
%-------------------------------------------------------------------------------
lgrid.main = uigridlayout(gui.Fig,[1 2],'Padding',[1 1 1 1]*spacing,'ColumnSpacing',spacing);
lgrid.main.ColumnWidth = {'7x','3x'};

% left panel: parameter table and plot
lgrid.left = uigridlayout(lgrid.main,[2 1],'Padding',[0 0 0 0]);
lgrid.left.RowHeight = {hPtop,'1x'};
lgrid.params = uigridlayout(lgrid.left,[2 1],'Padding',[0 0 0 0]);
lgrid.params.RowHeight = {hElement,'1x'};
lgrid.plot = uigridlayout(lgrid.left,[1 1]);

% right panel: settings
lgrid.right = uigridlayout(lgrid.main,[4 1],'Padding',[0 0 0 0],'RowSpacing',0.75*spacing);
lgrid.right.RowHeight = {hPtop,'1x','1x','1x'};

% Axes
%-------------------------------------------------------------------------------
% Data display
if numel(data)==1
  gui.datapanel = uigridlayout(lgrid.plot,[1 1],'Padding',[0 0 0 0]);
else
  lgrid.plot.Padding = [10 10 10 0];
  gui.datapanel = uigridlayout(lgrid.plot,[2 1],'Padding',[0 0 0 0]);
  gui.datapanel.RowHeight = {5.2*0.9*hElement, '1x'};

  multifitbox = uigridlayout(gui.datapanel,[1 2],'Padding',[0 0 0 0],'ColumnSpacing',20);
  multifitbox.ColumnWidth = {'1x','0.1x'};
  
  % Multifit table
  N = esfitdata.nDataSets;
  tablebox = uigridlayout(multifitbox,[6 2+N],'Padding',[0 0 0 0],'RowSpacing',0,'ColumnSpacing',0);
  tablebox.RowHeight = {0.9*hElement,0.9*hElement,0.9*hElement,0.9*hElement,'1x'};
  tmp = (0.7*sz(1)-4*spacing)/(5*hElement);
  if N>tmp
    tablebox.ColumnWidth = [repmat({'1x'},1,N+1) {'0x'}];
  else
    tablebox.ColumnWidth = [repmat({'1x'},1,N+1) {sprintf('%1.0fx',floor(tmp-N))}];
  end
  
  l(1) = uilabel('Parent',tablebox,...
                 'Text','Weight',...
                 'BackgroundColor',get(gui.Fig,'Color'),...
                 'FontWeight','bold',...
                 'HorizontalAl','left','VerticalAl','center');
  l(1).Layout.Column = 1;
  l(1).Layout.Row = 2;
  l(2) = uilabel('Parent',tablebox,...
                 'Text','Baseline',...
                 'BackgroundColor',get(gui.Fig,'Color'),...
                 'FontWeight','bold',...
                 'HorizontalAl','left','VerticalAl','center');
  l(2).Layout.Column = 1;
  l(2).Layout.Row = 3;  
  l(3) = uilabel('Parent',tablebox,...
                 'Text','Data points',...
                 'BackgroundColor',get(gui.Fig,'Color'),...
                 'FontWeight','bold',...
                 'HorizontalAl','left','VerticalAl','center');
  l(3).Layout.Column = 1;
  l(3).Layout.Row = 4;  
  l(4) = uilabel('Parent',tablebox,...
                 'Text','Scale',...
                 'BackgroundColor',get(gui.Fig,'Color'),...
                 'FontWeight','bold',...
                 'HorizontalAl','left','VerticalAl','center');
  l(4).Layout.Column = 1;
  l(4).Layout.Row = 5;  

  for i = 1:N

    % Dataset label
    c(i) = uilabel('Parent',tablebox,'Text',num2str(i),...
                   'FontWeight','bold','HorizontalAl','center','VerticalAl','center');
    c(i).Layout.Column = i+1;
    c(i).Layout.Row = 1;

    % Weight
    gui.multifitWeight(i) = uieditfield('text','Parent',tablebox,...
                                        'Value',num2str(esfitdata.Opts.weight(i)),...
                                        'Tooltip','dataset weight for global fitting (accepts numeric values and expressions, e.g. 1/240)',...
                                        'HorizontalAlignment','right',...
                                        'FontSize',fontsizetbl,...
                                        'ValueChangedFcn',@(src,evt) changemultifitWeight(src,evt));
    gui.multifitWeight(i).Layout.Column = i+1;
    gui.multifitWeight(i).Layout.Row = 2;
    
    % Baseline
    gui.multifitBaseline(i) = uidropdown('Parent',tablebox,...
                                  'Items',esfitdata.BaseLineStrings,...
                                  'ItemsData',esfitdata.BaseLineSettings,...
                                  'Value',esfitdata.BaseLine(i),...
                                  'FontSize',fontsizetbl,...
                                  'BackgroundColor','w',...
                                  'ValueChangedFcn',@(src,evt) changeBaseLineSelection(src,evt));
    
    gui.multifitBaseline(i).Layout.Column = i+1;
    gui.multifitBaseline(i).Layout.Row = 3;

    % Number of points
    gui.multifitDatasetpts(i) = uieditfield('numeric','Parent',tablebox,'Editable','off','FontSize',fontsizetbl);
    gui.multifitDatasetpts(i).Layout.Column = i+1;
    gui.multifitDatasetpts(i).Layout.Row = 4;
    gui.multifitDatasetpts(i).Value = numel(data{i});

    % Scale
    gui.multifitScale(i) = uieditfield('numeric','Parent',tablebox,'Editable','off','FontSize',fontsizetbl);
    gui.multifitScale(i).Layout.Column = i+1;
    gui.multifitScale(i).Layout.Row = 5;
    
  end

  % Tile vs tab display toggle buttons
  buttonbox = uigridlayout(multifitbox,[4 1],'Padding',[0 0 0 0],'RowSpacing',0);
  buttonbox.RowHeight = {hElement,hElement,hElement,'1x'};
  gui.TilesButton = uibutton('state','Parent',buttonbox,...
                             'Text','Tiles','Enable','on',...
                             'ValueChangedFcn',@(src,evt) setDataDisplayType(src,'Tiles'),...
                             'Tooltip','Show datasets in different tiles');
  gui.TilesButton.Layout.Row = 2;
  gui.TabsButton = uibutton('state','Parent',buttonbox,...
                             'Text','Tabs','Enable','on',...
                             'ValueChangedFcn',@(src,evt) setDataDisplayType(src,'Tabs'),...
                             'Tooltip','Show datasets in different tabs');
  gui.TabsButton.Layout.Row = 3;

  if esfitdata.nDataSets<=8
    gui.TilesButton.Value = 1;
    gui.TabsButton.Value = 0;
  else
    gui.TilesButton.Value = 0;
    gui.TabsButton.Value = 1;
  end

end

setupDataDisplay('create');

showmaskedregions();

% Fit parameter display
%-------------------------------------------------------------------------------
lgrid.parbanner = uigridlayout(lgrid.params,[1 6],'Padding',[0 0 0 0],'ColumnSpacing',20);
lgrid.parbanner.ColumnWidth = {'2x','4.25x','0.1x','4.25x','0.1x','3x'};
lgrid.parbanner.Layout.Row = 1;
uilabel('Parent',lgrid.parbanner,...
        'Text','Parameters',...
        'BackgroundColor',get(gui.Fig,'Color'),...
        'FontWeight','bold',...
        'HorizontalAl','left','VerticalAl','center');

startpointbox = uigridlayout(lgrid.parbanner,[1 4],'Padding',[0 0 0 0],'ColumnSpacing',0);
startpointbox.ColumnWidth = {'1.25x','1x','1x','1x'};
uilabel('Parent',startpointbox,...
    'BackgroundColor',get(gui.Fig,'Color'),...
    'FontWeight','bold','Text','Start point:',...
    'HorizontalAl','left','VerticalAl','center');
gui.selectStartPointButtonCenter = uibutton('Parent',startpointbox,...
    'Text','center','Enable','on','ButtonPushedFcn',@(src,evt) setStartPoint('center'),...
    'Tooltip','Set start values to center of range');
gui.selectStartPointButtonRandom = uibutton('Parent',startpointbox,...
    'Text','random','Enable','on','ButtonPushedFcn',@(src,evt) setStartPoint('random'),...
    'Tooltip','Set random start values');
gui.selectStartPointButtonBest = uibutton('Parent',startpointbox,...
    'Text','best','Enable','on','ButtonPushedFcn',@(src,evt) setStartPoint('best'),...
    'Tooltip','Set start values to current best fit');
  
selectbox = uigridlayout(lgrid.parbanner,[1 4],'Padding',[0 0 0 0],'ColumnSpacing',0);
selectbox.ColumnWidth = startpointbox.ColumnWidth;
selectbox.Layout.Column = 4;
uilabel('Parent',selectbox,...
    'BackgroundColor',get(gui.Fig,'Color'),...
    'FontWeight','bold','Text','Selection:',...
    'HorizontalAl','left','VerticalAl','center');
gui.selectInvButton = uibutton('Parent',selectbox,...
    'Text','invert','Enable','on','ButtonPushedFcn',@(src,evt) selectButtonCallback('invert'),...
    'Tooltip','Invert selection of parameters');
gui.selectAllButton = uibutton('Parent',selectbox,...
    'Text','all','Enable','on','ButtonPushedFcn',@(src,evt) selectButtonCallback('all'),...
    'Tooltip','Select all parameters');
gui.selectNoneButton = uibutton('Parent',selectbox,...
    'Text','none','Enable','on','ButtonPushedFcn',@(src,evt) selectButtonCallback('none'),...
    'Tooltip','Unselect all parameters');

fitparnumbox = uigridlayout(lgrid.parbanner,[1 2],'Padding',[0 0 0 0],'ColumnSpacing',0);
fitparnumbox.ColumnWidth = {'2x','1x'};
fitparnumbox.Layout.Column = 6;
uilabel('Parent',fitparnumbox,...
    'BackgroundColor',get(gui.Fig,'Color'),...
    'FontWeight','bold','Text','Active parameters:',...
    'HorizontalAl','left','VerticalAl','center');
gui.nParamsField = uieditfield('numeric','Parent',fitparnumbox,...
    'Editable','off','FontSize',fontsizetbl,...
    'Value',numel(esfitdata.pinfo),...
    'Tooltip','Number of selected active parameters');

% Parameter table
columnname = {'','','Name','start','lower','upper','current','best','ci95 lower','ci95 upper'};
columnformat = {'char','logical','char','char','char','char','char','char','char','char'};
colEditable = [false true false true true true false false false false];
data = cell(numel(esfitdata.pinfo),10);
for p = 1:numel(esfitdata.pinfo)
  data{p,1} = num2str(p);
  data{p,2} = true;
  data{p,3} = char(esfitdata.pinfo(p).Name);
  data{p,4} = sprintf('%0.6g',esfitdata.p_start(p));
  data{p,5} = sprintf('%0.6g',esfitdata.pvec_lb(p));
  data{p,6} = sprintf('%0.6g',esfitdata.pvec_ub(p));
  data{p,7} = '-';
  data{p,8} = '-';
  data{p,9} = '-';
  data{p,10} = '-';
end
columnwidths = {hElement,hElement,'auto','auto','auto','auto','auto','auto','auto','auto'};
gui.ParameterTable = uitable('Parent',lgrid.params,...
                             'ColumnFormat',columnformat,'ColumnName',columnname,'RowName',[],...
                             'ColumnEditable',colEditable,'FontSize',fontsizetbl,...
                             'CellEditCallback',@paramsTableCellEditCallback,...
                             'ColumnWidth',columnwidths,...
                             'Data',data,...
                             'UserData',colEditable);
gui.ParameterTable.Layout.Row = 2;
if ~verLessThan('Matlab','9.13')
  s = uistyle('Interpreter','html');
  addStyle(gui.ParameterTable,s);
end

% FitOption selection
%-------------------------------------------------------------------------------
lgrid.fitcontrol = uigridlayout(lgrid.right,[1 2],'Padding',[0 0 0 0],'ColumnSpacing',spacing);
lgrid.fitcontrol.ColumnWidth = {'1.25x','1x'};

fitoptbox0 = uigridlayout(lgrid.fitcontrol,[4 1],'Padding',[0 0 0 0],'ColumnSpacing',0,'RowSpacing',5);
fitoptbox0.RowHeight = {hElement,'4x','3x'};

fitoptbox1 = uigridlayout(fitoptbox0,[1 2],'Padding',[0 0 0 0],'ColumnSpacing',0);
fitoptbox1.ColumnWidth = {'0.6x','0.8x',hElement};
uilabel('Parent',fitoptbox1,...
    'Text','Function',...
    'Tooltip','Name of simulation function',...
    'FontWeight','bold',...
    'HorizontalAl','left','VerticalAl','center',...
    'BackgroundColor',get(gui.Fig,'Color'));
uilabel('Parent',fitoptbox1,...
    'Text',esfitdata.fcnName,...
    'FontColor','b',...
    'Tooltip',{esfitdata.fcnName,sprintf('using output no. %d of %d',esfitdata.OutArgument,esfitdata.nOutArguments)},...
    'HorizontalAl','left','VerticalAl','center',...
    'BackgroundColor',get(gui.Fig,'Color'));
gui.AlgorithmSettingsButton = uibutton('Parent',fitoptbox1,...
                                       'Text','','Tooltip','Set fitting algorithm options',...
                                       'Icon','private/settingsicon.png','IconAlignment','center',...
                                       'Enable','on');

fitoptbox2 = uigridlayout(fitoptbox0,[4 2],'Padding',[0 0 0 0],'ColumnSpacing',0,'RowSpacing',5);
fitoptbox2.ColumnWidth = {'0.6x','1x'};
uilabel('Parent',fitoptbox2,...
    'Text','Algorithm',...
    'FontWeight','bold',...
    'HorizontalAlign','left','VerticalAlign','center',...
    'BackgroundColor',get(gui.Fig,'Color'));
gui.AlgorithmMenu = uidropdown('Parent',fitoptbox2,...
                               'Items',esfitdata.AlgorithmNames,...
                               'ItemsData',1:numel(esfitdata.AlgorithmNames),...
                               'Value',Opt.AlgorithmID,...
                               'BackgroundColor','w',...
                               'ValueChangedFcn',@(src,evt) selectAlgorithm(src),...
                               'Tooltip','Fitting algorithm');

uilabel('Parent',fitoptbox2,...
    'Text','Target',...
    'FontWeight','bold',...
    'HorizontalAlign','left','VerticalAlign','center',...
    'BackgroundColor',get(gui.Fig,'Color'));
gui.TargetMenu = uidropdown('Parent',fitoptbox2,...
                            'Items',esfitdata.TargetNames,...
                            'ItemsData',1:numel(esfitdata.TargetNames),...
                            'Value',Opt.TargetID,...
                            'BackgroundColor','w',...
                            'Tooltip','Target function');

uilabel('Parent',fitoptbox2,...
    'Text','AutoScale',...
    'FontWeight','bold',...
    'HorizontalAlign','left','VerticalAlign','center',...
    'BackgroundColor',get(gui.Fig,'Color'));
gui.AutoScaleMenu = uidropdown('Parent',fitoptbox2,...
                               'Items',esfitdata.AutoScaleStrings,...
                               'ItemsData',esfitdata.AutoScaleSettings,...
                               'Value',esfitdata.Opts.AutoScaleID,...
                               'BackgroundColor','w',...
                               'Tooltip','Autoscaling');

uilabel('Parent',fitoptbox2,...
    'Text','BaseLine',...
    'FontWeight','bold',...
    'HorizontalAlign','left','VerticalAlign','center',...
    'BackgroundColor',get(gui.Fig,'Color'));
gui.BaseLineMenu = uidropdown('Parent',fitoptbox2,...
                              'Items',esfitdata.BaseLineStrings,...
                              'ItemsData',esfitdata.BaseLineSettings,...
                              'BackgroundColor','w');
if esfitdata.nDataSets==1
  gui.BaseLineMenu.Value = esfitdata.BaseLine;
  gui.BaseLineMenu.Tooltip = 'Baseline fitting';
else
  set(gui.BaseLineMenu,'Value',esfitdata.BaseLine(1));
  set(gui.BaseLineMenu,'Tooltip','Baseline fitting (global setting applied to all data sets)');
  set(gui.BaseLineMenu,'ValueChangedFcn',@(src,evt) changeBaseLineSelection(src,evt));
end

fitoptbox3 = uigridlayout(fitoptbox0,[3 2],'Padding',[0 0 0 0],'ColumnSpacing',5,'RowSpacing',5);
gui.MaskCheckbox = uicheckbox('Parent',fitoptbox3,...
           'Text','use mask','Tooltip','Use mask with excluded regions',...
           'Value',1);
gui.clearMaskButton = uibutton('Parent',fitoptbox3,...
         'Text','clear mask','Tooltip','Clear mask',...
         'Enable','on',...
         'ButtonPushedFcn',@clearMaskCallback);

gui.ErrorCalcCheckbox = uicheckbox('Parent',fitoptbox3,...
           'Text','calculate uncertainties','Tooltip','Calculate uncertainties for fitting parameters after convergence',...
           'Value',esfitdata.Opts.CalculateUncertainties,'ValueChangedFcn','global gui esfitdata; esfitdata.Opts.CalculateUncertainties = gui.ErrorCalcCheckbox.Value;');
gui.ErrorCalcCheckbox.Layout.Column = [1 2];

gui.BaselineCheckbox = uicheckbox('Parent',fitoptbox3,...
           'Text','plot baseline','Tooltip','Plot baseline on top of data',...
           'Value',0,'ValueChangedFcn',@showbaseline);
gui.ResidualCheckbox = uicheckbox('Parent',fitoptbox3,...
           'Text','plot residuals','Tooltip','Plot residuals',...
           'Value',0,'ValueChangedFcn',@showresiduals);

% Start/Stop buttons
%--------------------------------------------------------------------------
fitctrlbox = uigridlayout(lgrid.fitcontrol,[4 1],'Padding',[0 0 0 0],'RowSpacing',0);
fitctrlbox.RowHeight = {'2x','1x','1x','1x'};
gui.StartButton = uibutton('Parent',fitctrlbox,...
                           'Text','Start fitting',...
                           'ButtonPushedFcn',@startButtonCallback,...
                           'Visible','on',...
                           'Tooltip','Start fitting');
gui.SaveButton = uibutton('Parent',fitctrlbox,...
                          'Text','Save parameter set',...
                          'ButtonPushedFcn',@saveFitsetCallback,...
                          'Enable','off',...
                          'Tooltip','Save latest fitting result');
gui.EvaluateButton = uibutton('Parent',fitctrlbox,...
                              'Text','Evaluate at start point',...
                              'ButtonPushedFcn',@evaluateCallback,...
                              'Enable','on',...
                              'Tooltip','Run simulation for current start parameters');
gui.ResetButton = uibutton('Parent',fitctrlbox,...
                           'Text','Reset',...
                           'ButtonPushedFcn',@resetCallback,...
                           'Enable','on',...
                           'Tooltip','Clear fit history');

% Fitset list
%-------------------------------------------------------------------------------
lgrid.fitsets = uigridlayout(lgrid.right,[2 1],'Padding',[0 0 0 0],'RowSpacing',5);
lgrid.fitsets.RowHeight = {hElement,'1x',hElement};

fitsettitlebox = uigridlayout(lgrid.fitsets,[1 2],'Padding',[0 0 0 0],'RowSpacing',5);
fitsettitlebox.ColumnWidth = {'1x',hElement};
uilabel('Parent',fitsettitlebox,...
    'BackgroundColor',get(gui.Fig,'Color'),...
    'FontWeight','bold','Text','Fit parameter sets',...
    'Tooltip','List of stored fit parameter sets',...
    'HorizontalAl','left','VerticalAl','top');
gui.autosave = uibutton('state','Parent',fitsettitlebox,...
         'Text','','Tooltip','Turn autosave on/off',...
         'Icon','private/saveicon.png','IconAlignment','center',...
         'Enable','on','Value',false);


% Fit parameter set table
columnname = {'ID','rmsd','algorithm','target','scale','baseline',''};
columnformat = {'char','char','char','char','char','char','char'};
colEditable = [false false false false false false false];
columnwidths = {1.6*hElement,'auto','auto','auto','auto','auto','auto'};
gui.FitSetTable = uitable('Parent',lgrid.fitsets,...
                          'ColumnFormat',columnformat,'ColumnName',columnname,'RowName',[],...
                          'ColumnEditable',colEditable,'FontSize',fontsizetbl,...
                          'ColumnWidth',columnwidths,...
                          'ColumnSortable',true,...
                          'Data',[],...
                          'SelectionType','row',...
                          'CellSelectionCallback',@setListCallback,...
                          'KeyPressFcn',@deleteSetListKeyPressFcn,...
                          'Enable','on');
fitsetbuttonbox = uigridlayout(lgrid.fitsets,[1 5],'Padding',[0 0 0 0],'ColumnSpacing',0);
fitsetbuttonbox.ColumnWidth = {'1x','1x','1x','1x','1x'};
gui.selectStartPointButtonSelected = uibutton('Parent',fitsetbuttonbox,...
                                              'Text','set as start',...
                                              'Tooltip','Set as start point for fitting','Enable','off',...
                                              'ButtonPushedFcn',@(src,evt) setStartPoint('selected'));
gui.selectStartPointButtonSelected.Layout.Column = 1;
gui.exportSetButton = uibutton('Parent',fitsetbuttonbox,...
    'Text','export',...
    'Tooltip','Export fit parameter set to workspace','Enable','off',...
    'ButtonPushedFcn',@exportSetButtonCallback);
gui.exportSetButton.Layout.Column = 2;
gui.saveSetButton = uibutton('Parent',fitsetbuttonbox,...
    'Text','save',...
    'Tooltip','Save fit parameter set as mat file','Enable','off',...
    'ButtonPushedFcn',@saveSetButtonCallback);
gui.saveSetButton.Layout.Column = 3;
gui.deleteSetButton = uibutton('Parent',fitsetbuttonbox,...
    'Text','delete',...
    'Tooltip','Delete fit parameter set','Enable','off',...
    'ButtonPushedFcn',@deleteSetButtonCallback);
gui.deleteSetButton.Layout.Column = 4;
gui.clearSetButton = uibutton('Parent',fitsetbuttonbox,...
    'Text','clear',...
    'Tooltip','Clear fit parameter set table','Enable','off',...
    'ButtonPushedFcn',@clearSetButtonCallback);
gui.clearSetButton.Layout.Column = 5;

% Iteration and rmsd history displays
%-------------------------------------------------------------------------------
lgrid.rmsdplot = uigridlayout(lgrid.right,[4 1],'Padding',[0 0 0 0],'RowSpacing',0);
lgrid.rmsdplot.RowHeight = {hElement,hElement,'1x',hElement};

uilabel('Parent',lgrid.rmsdplot,...
    'BackgroundColor',get(gui.Fig,'Color'),...
    'FontWeight','bold','Text','RMSD history',...
    'HorizontalAl','left','VerticalAl','center');

rmsdinfobox = uigridlayout(lgrid.rmsdplot,[1 2],'Padding',[0 0 0 0]);
rmsdinfobox.ColumnWidth = {'1x',hElement};
gui.RmsText = uilabel('Parent',rmsdinfobox,...
        'Text',' RMSD: -',...
        'FontColor',[0 0.6 0],'Tooltip','Current best RMSD',...
        'HorizontalAl','left','VerticalAl','center');

gui.RmsLogPlot = uibutton('state','Parent',rmsdinfobox,...
         'Text','log','FontSize',fontsizetbl-1,'Tooltip','Set log scale on/off',...
         'Value',0,'ValueChangedFcn',@updatermsdplot);

hAx = uiaxes('Parent',lgrid.rmsdplot,'Units','pixels','Layer','top');
gui.rmsdline = line(hAx,1,NaN,'Marker','.','LineStyle','none','DisplayName','rmsdline');
set(gui.rmsdline,'MarkerSize',5,'Color',[0.2 0.2 0.8]);
set(hAx,'box','on','YScale','lin','XTick',[],'YAxisLoc','right','Layer','top','YGrid','on');

gui.logLine = uilabel('Parent',lgrid.rmsdplot,...
          'Text','','WordWrap','on',...
          'FontSize',fontsizetbl,...
          'Tooltip','Information from fitting algorithm',...
          'HorizontalAl','left','VerticalAl','center');

% Error log panel
%-------------------------------------------------------------------------------
lgrid.log = uigridlayout(lgrid.right,[2 1],'Padding',[0 0 0 0],'RowSpacing',5);
lgrid.log.RowHeight = {hElement,'1x'};

logtitlebox = uigridlayout(lgrid.log,[1 2],'Padding',[0 0 0 0],'RowSpacing',5);
logtitlebox.ColumnWidth = {'1x',hElement};
uilabel('Parent',logtitlebox,...
    'BackgroundColor',get(gui.Fig,'Color'),...
    'FontWeight','bold','Text','Log',...
    'Tooltip','Fitting information and error log',...
    'HorizontalAl','left','VerticalAl','top');
uibutton('Parent',logtitlebox,...
         'Text','','Tooltip','Copy to clipboard',...
         'Icon','private/copyicon.png','IconAlignment','center',...
         'Enable','on',...
         'ButtonPushedFcn',@copyLog);

gui.LogBox = uitextarea('Parent',lgrid.log,...
                        'Value',{''},'Tooltip','',...
                        'WordWrap','on',...
                        'FontSize',fontsizetbl,...
                        'Editable','off',...
                        'BackgroundColor',[1 1 1]);

copymenu = uicontextmenu(gui.Fig);
uimenu(copymenu,'Text','Copy to clipboard','Callback',@copyLog);
gui.LogBox.ContextMenu = copymenu;
drawnow limitrate

% Algorithm settings pop-up figure
%-------------------------------------------------------------------------------
gui.Popup = uifigure('Tag','algorithmpopup','WindowStyle','normal','Visible','off');

sz = [400 320]; % popup size
sz = sz*scalefact;
xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the screen vertically
set(gui.Popup,'position',[xpos, ypos, sz(1), sz(2)],'units','pixels');
set(gui.Popup,'MenuBar','none');
% set(gui.Popup,'Resize','off');
set(gui.Popup,'Name','esfit - Algorithm settings','NumberTitle','off');
set(gui.Popup,'CloseRequestFcn',@openAlgorithmSettings);
  
% Set up drop down menu and tabs
%-------------------------------------------------------------------------------
lgrid.main = uigridlayout(gui.Popup,[3 1],'Padding',[1 1 1 1]*spacing,'RowSpacing',5);
lgrid.main.RowHeight = {hElement,hElement,'1x'};

dropdownbox = uigridlayout(lgrid.main,[1 4],'Padding',[0 0 0 0],'ColumnSpacing',0);
dropdownbox.ColumnWidth = {'0.6x','1x','0.2x','0.6x'};
uilabel('Parent',dropdownbox,...
    'Text','Algorithm',...
    'FontWeight','bold',...
    'HorizontalAlign','left','VerticalAlign','center',...
    'BackgroundColor',get(gui.Popup,'Color'));
gui.AlgorithmMenuPopup = uidropdown('Parent',dropdownbox,...
    'Items',esfitdata.AlgorithmNames,...
    'ItemsData',1:numel(esfitdata.AlgorithmNames),...
    'Value',Opt.AlgorithmID,...
    'BackgroundColor','w',...
    'ValueChangedFcn',@(src,evt) selectAlgorithm(src),...
    'Tooltip','Fitting algorithm');
gui.setAlgorithmDefaults = uibutton('Parent',dropdownbox,...
         'Text','Set to defaults','Enable','on','ButtonPushedFcn',@(src,evt) updateAlgorithmDefaults,...
         'Tooltip','Reset all parameters to the default');
gui.setAlgorithmDefaults.Layout.Column = 4;

nparamsbox = uigridlayout(lgrid.main,[1 4],'Padding',[0 0 0 0],'ColumnSpacing',0);
nparamsbox.ColumnWidth = {'1.2x','0.4x','0.4x','1x'};
l = uilabel('Parent',nparamsbox,...
    'Text','N',...
    'FontWeight','bold',...
    'HorizontalAlign','center','VerticalAlign','center',...
    'BackgroundColor',get(gui.Popup,'Color'));
l.Layout.Column = 2;
gui.nParamsFieldPopup = uieditfield('numeric','Parent',nparamsbox,...
        'ToolTip','number of fitting parameters',...
        'Value',esfitdata.nParameters,...
        'Editable','off');

gui.AlgorithmTabs = uitabgroup(lgrid.main,'TabLocation','left');
for i = 1:(numel(esfitdata.AlgorithmNames)-1)
  t(i) = uitab(gui.AlgorithmTabs,'Title',esfitdata.AlgorithmAbbrev{i},'ToolTip',esfitdata.AlgorithmNames{i});
  [FitOpt,info] = esfit_algdefaults(esfitdata.AlgorithmNames{i});
  esfitdata.AlgorithmDefaults{i} = FitOpt;
  setting = fieldnames(info);
  tabgrid = uigridlayout(t(i),[7 2],'Padding',[1 1 1 1]*spacing/2,'RowSpacing',5,'ColumnSpacing',spacing);
  tabgrid.RowHeight = [num2cell(repmat(hElement,1,numel(setting))),{'1x'}];
  tabgrid.ColumnWidth = {'0.8x','1x'};
  clear hsetting
  for j = 1:numel(setting)
    uilabel('Parent',tabgrid,...
        'Text',setting{j},...
        'FontWeight','bold',...
        'HorizontalAlign','left','VerticalAlign','center',...
        'BackgroundColor',get(gui.Popup,'Color'));
    defaultpar = FitOpt.(setting{j});
    if islogical(defaultpar)
      hsetting(j) = uidropdown('Parent',tabgrid,...
          'Items',{'true','false'},...
          'ItemsData',[1 0],...
          'ToolTip',info.(setting{j}),...
          'BackgroundColor','w',...
          'UserData',{setting{j},'logical'},...
          'ValueChangedFcn',@(src,evt) changeAlgorithmSetting(src,evt));
    elseif isnumeric(defaultpar) && numel(defaultpar)==1
      hsetting(j) = uieditfield('Parent',tabgrid,...
        'ToolTip',info.(setting{j}),...
        'UserData',{setting{j},'num'},...
        'ValueChangedFcn',@(src,evt) changeAlgorithmSetting(src,evt));
    else
      hsetting(j) = uieditfield('Parent',tabgrid,...
        'ToolTip',info.(setting{j}),...
        'UserData',{setting{j},'eval'},...
        'ValueChangedFcn',@(src,evt) changeAlgorithmSetting(src,evt));
    end
  end
  gui.AlgorithmTabs.UserData{i} = hsetting;
end
updateAlgorithmDefaults()

% Set callback for settings button to open popup menu
settingsbutton = gui.AlgorithmSettingsButton;
set(settingsbutton,'ButtonPushedFcn',@openAlgorithmSettings);

set(gui.Fig,'AutoResizeChildren',false)
set(gui.Popup,'AutoResizeChildren',false)
set(gui.Fig,'SizeChangedFcn',@(src,evt) resizeGUI('main'))
set(gui.Popup,'SizeChangedFcn',@(src,evt) resizeGUI('popup'))
set(gui.Fig,'Visible','on')
set(gui.Fig,'NextPlot','new');
set(gui.Popup,'NextPlot','new');

end
%===============================================================================

%===============================================================================
function setDataDisplayType(src,type)

global gui

statuschanged = false;

switch type
  case 'Tiles'
    if src.Value==1
      statuschanged = true;
    end
    gui.TilesButton.Value = 1;
    gui.TabsButton.Value = 0;
    %   gui.TabsButton.Value = 0;
    % else
    %   gui.TabsButton.Value = 1;
    % end
  case 'Tabs'
    if src.Value==1
      statuschanged = true;
    end
    gui.TabsButton.Value = 1;
    gui.TilesButton.Value = 0;
    % if src.Value==1
      % gui.TilesButton.Value = 0;
    % else
      % gui.TilesButton.Value = 1;
    % end
end

if statuschanged
  setupDataDisplay('switch')
end

end
%===============================================================================

%===============================================================================
function setupDataDisplay(evt)

global esfitdata gui
Opt = esfitdata.Opts;

N = esfitdata.nDataSets;

if N==1
  if strcmp(evt,'create')
    gui.dataaxes = uiaxes('Parent',gui.datapanel,'Layer','top','Tag','1');
  end
else

  if strcmp(evt,'switch')

    % Get current visibility of plot components
    fieldname = {'currsimdata','bestsimdata','residualdata','residualzero','baselinedata'};
    showdata = zeros(1,numel(fieldname));
    for j = 1:numel(fieldname)
      if gui.(fieldname{j})(1).Visible
        showdata(j) = 1;
      end
    end

    delete(gui.datapanel.Children(2))

  end

  if gui.TilesButton.Value
    % Axes in tiles
    nx = floor(sqrt(N));
    ny = ceil(N/nx);

    databox = uigridlayout(gui.datapanel,[nx ny],'Padding',[0 0 0 0],'Tag','data');
    set(databox,'Visible','off');
    for i = 1:N
      gui.dataaxes(i) = uiaxes('Parent',databox,'Layer','top','Tag',num2str(i));
    end
    drawnow
    set(databox,'Visible','on');
  elseif gui.TabsButton.Value
    % Axes in tabs
    axestabs = uitabgroup(gui.datapanel,'TabLocation','top','Tag','data');
    set(axestabs,'Visible','off');
    for i = 1:N
      gui.datatab(i) = uitab(axestabs,'Title',num2str(i));
      tabgrid = uigridlayout(gui.datatab(i),[1 1]);
      gui.dataaxes(i) = uiaxes('Parent',tabgrid,'Layer','top','Tag',num2str(i));
    end
    drawnow
    set(axestabs,'Visible','on');
  end
end

for i = 1:N

  expdata = esfitdata.data(esfitdata.idx{i});

  NaNdata = NaN(1,numel(expdata));
  
  mask = esfitdata.Opts.Mask(esfitdata.idx{i});
  maxy = max(expdata(mask));
  miny = min(expdata(mask));
  YLimits = [miny maxy] + [-1 1]*Opt.PlotStretchFactor*(maxy-miny);
  x = esfitdata.Opts.x(esfitdata.idx{i});
  minx = min(x);
  maxx = max(x);

  gui.baselinedata(i) = line(gui.dataaxes(i),x,NaNdata,'Color',[1 1 1]*0.65,'LineWidth',1,'Visible','off','DisplayName','baseline');
  gui.expdata(i) = line(gui.dataaxes(i),x,expdata,'Color','k','Marker','.','LineStyle','none','DisplayName','expdata');
  gui.currsimdata(i) = line(gui.dataaxes(i),x,NaNdata,'Color','r','DisplayName','currsim');
  gui.bestsimdata(i) = line(gui.dataaxes(i),x,NaNdata,'Color',[0 0.6 0],'DisplayName','bestsim');
  gui.dataaxes(i).XLim = [minx maxx];
  gui.dataaxes(i).YLim = YLimits;
  gui.dataaxes(i).ButtonDownFcn = @(src,evt) axesButtonDownFcn(src);
  grid(gui.dataaxes(i),'on');
  box(gui.dataaxes(i),'on')
  if ~gui.showxaxis % no x-axis defined
    set(gui.dataaxes(i),'XTickLabel',{})
  end
  set(gui.dataaxes(i),'YTickLabel',{})

  if isfield(gui,'TilesButton') && gui.TilesButton.Value
    gui.tileID(i) = text(gui.dataaxes(i),minx+0.05*(maxx-minx),YLimits(2)-0.05*diff(YLimits),sprintf('%i',i),'FontWeight','bold');
  end

  gui.residualzero(i) = yline(gui.dataaxes(i),miny-abs((maxy-miny)/3),'DisplayName','residualzero');
  gui.residualdata(i) = line(gui.dataaxes(i),x,NaNdata,'Color',[1 1 1]*0.65,'Marker','.','LineStyle','none','DisplayName','residuals');
  set(gui.residualzero(i),'Visible','off')
  set(gui.residualdata(i),'Visible','off')
end

if strcmp(evt,'switch')

  for i = 1:N

    idx = esfitdata.idx{i};

    % Get relevant quantities
    x = esfitdata.Opts.x(idx);
    expdata = esfitdata.data(idx);

    if isfield(esfitdata,'curr') && isfield(esfitdata.curr,'sim')
      currsim = real(esfitdata.curr.sim(idx));
      currbaseline = real(esfitdata.curr.baseline(idx));
    else
      currsim = NaN*ones(size(expdata));
      currbaseline = NaN*ones(size(expdata));
    end
    if isfield(esfitdata,'best') && isfield(esfitdata.best,'fit')
      bestsim = real(esfitdata.best.fit(idx));
    else
      bestsim = NaN*ones(size(expdata));
    end

    % Update plotted data
    set(gui.expdata(i),'XData',x,'YData',expdata);
    set(gui.bestsimdata(i),'XData',x,'YData',bestsim);
    set(gui.currsimdata(i),'XData',x,'YData',currsim);
    set(gui.baselinedata(i),'XData',x,'YData',currbaseline);
    residualzero = get(gui.residualzero(i),'Value');
    residualdata = (expdata(:)-currsim(:))+residualzero;
    set(gui.residualdata(i),'XData',x,'YData',residualdata);

    % Update visibility
    for j = 1:numel(fieldname)
      if showdata(j)==1
        set(gui.(fieldname{j})(i),'Visible','on');
      else
        set(gui.(fieldname{j})(i),'Visible','off');
      end
    end

  end

  showmaskedregions()
  updateaxislimits()

end

end
%===============================================================================

%===============================================================================
function evaluateCallback(~,~)
% Evaluate for selected parameters
global esfitdata gui

esfitdata.modelErrorHandler = @(ME) GUIErrorHandler(ME);

p_eval = esfitdata.p_start;
active = ~esfitdata.fixedParams;
p_eval = p_eval(active);
esfitdata.Opts.AutoScaleID = get(gui.AutoScaleMenu,'Value');
if esfitdata.nDataSets==1
  esfitdata.Opts.BaseLine = get(gui.BaseLineMenu,'Value');
else
  weight = zeros(1,sum(esfitdata.datasize));
  for i = 1:esfitdata.nDataSets
    esfitdata.Opts.BaseLine(i) = get(gui.multifitBaseline(i),'Value');
    esfitdata.Opts.weight(i) = str2double(get(gui.multifitWeight(i),'Value'));
    weight(esfitdata.idx{i}) = esfitdata.Opts.weight(i);
  end
  esfitdata.weight = weight;
end
esfitdata.Opts.useMask = get(gui.MaskCheckbox,'Value')==1;
Opt = esfitdata.Opts;
Opt.track = false;

try
  [~,rmsd] = residuals_(p_eval,Opt,1);
catch ME
  if isfield(esfitdata,'modelEvalError') && esfitdata.modelEvalError
    return
  elseif contains(ME.stack(1).name,'esfit')
    GUIErrorHandler(ME);
    return;
  else
    error(ME.message)
  end
end

for i = 1:esfitdata.nDataSets

  idx = esfitdata.idx{i};

  % Get current spectrum
  expdata = esfitdata.data(idx);
  currsim = real(esfitdata.curr.sim(idx));
  currbaseline = real(esfitdata.curr.baseline(idx));

  % Update plotted data
  x = esfitdata.Opts.x(idx);
  set(gui.expdata(i),'XData',x,'YData',expdata);
  set(gui.currsimdata(i),'XData',x,'YData',currsim);
  set(gui.baselinedata(i),'XData',x,'YData',currbaseline);
  residualzero = get(gui.residualzero(i),'Value');
  residualdata = (expdata(:)-currsim(:))+residualzero;
  set(gui.residualdata(i),'XData',x,'YData',residualdata);

  if esfitdata.nDataSets>1
    gui.multifitScale(i).Value = esfitdata.curr.scale(i);
  end

end

% Readjust vertical range
showmaskedregions()
updateaxislimits()
drawnow

% Update column with best values if current parameter set is new best
str = sprintf('Current RMSD: %g\n',esfitdata.curr.rmsd);
set(gui.RmsText,'Text',str,'FontColor',[1 0 0]);

end
%===============================================================================

%===============================================================================
function startButtonCallback(~,~)

global esfitdata gui

switch get(gui.StartButton,'Text')
  case 'Start fitting'
    esfitdata.UserCommand = 0;
  case 'Stop fitting'
    esfitdata.UserCommand = 1;
    return
end

% Check that there are selected fitting parameters
Data = gui.ParameterTable.Data;
if ~any([Data{:,2}])
  % Skip fitting and evaluate fit function with selected parameters
  msg = 'No active parameters; skipping optimization';
  esfitdata.UserCommand = 1;
  updateLogBox(msg);
  evaluateCallback()
  return;
end

% Update GUI
%-------------------------------------------------------------------------------

% Hide Start button, show Stop button
set(gui.StartButton,'Text','Stop fitting');
set(gui.StartButton,'Tooltip','Stop fitting');
set(gui.SaveButton,'Enable','off');

% Disable other buttons
set(gui.EvaluateButton,'Enable','off');
set(gui.ResetButton,'Enable','off');
set(gui.AlgorithmSettingsButton,'Enable','off');

% Disable listboxes
set(gui.AlgorithmMenu,'Enable','off');
set(gui.TargetMenu,'Enable','off');
set(gui.AutoScaleMenu,'Enable','off');
set(gui.BaseLineMenu,'Enable','off');
if esfitdata.nDataSets>1
  set(gui.multifitBaseline,'Enable','off');
  set(gui.multifitWeight,'Enable','off');
end

% Disable parameter table
set(gui.selectAllButton,'Enable','off');
set(gui.selectNoneButton,'Enable','off');
set(gui.selectInvButton,'Enable','off');
set(gui.selectStartPointButtonCenter,'Enable','off');
set(gui.selectStartPointButtonRandom,'Enable','off');
set(gui.selectStartPointButtonBest,'Enable','off');
colEditable = get(gui.ParameterTable,'UserData');
set(gui.ParameterTable,'ColumnEditable',false(size(colEditable)));
set(gui.ParameterTable,'CellEditCallback',[]);

% Disable algorithm settings window
set(gui.AlgorithmMenuPopup,'Enable','off');
set(gui.setAlgorithmDefaults,'Enable','off');
nTabs = numel(gui.AlgorithmTabs.UserData);
for i = 1:nTabs
  hsetting = gui.AlgorithmTabs.UserData{i};
  for j = 1:numel(hsetting)
    set(hsetting(j),'Enable','off');
  end
end

% Disable multifit settings
if esfitdata.nDataSets>1
  % set(gui.TilesButton,'Enable','off');
  % set(gui.TabsButton,'Enable','off');
  set(gui.multifitWeight,'Enable','off');
  set(gui.multifitBaseline,'Enable','off');
end

% Remove displayed best fit and uncertainties
Data = gui.ParameterTable.Data;
for p = 1:size(Data,1)
  Data{p,7} = '-';
  Data{p,8} = '-';
  Data{p,9} = '-';
  Data{p,10} = '-';
end
set(gui.ParameterTable,'Data',Data);

% Get fixed parameters
for p = 1:size(Data,1)
  esfitdata.fixedParams(p) = Data{p,2}==0;
end

% Disable fitset list controls
set(gui.selectStartPointButtonSelected,'Enable','off');
set(gui.exportSetButton,'Enable','off');
set(gui.saveSetButton,'Enable','off');
set(gui.deleteSetButton,'Enable','off');
set(gui.clearSetButton,'Enable','off');

% Disable mask tools
for i = 1:numel(gui.dataaxes)
  gui.dataaxes(i).ButtonDownFcn = [];
end
set(gui.clearMaskButton,'Enable','off');
set(gui.MaskCheckbox,'Enable','off');

% Pull settings from UI
%-------------------------------------------------------------------------------
% Determine selected method, target, autoscaling, start point
esfitdata.Opts.AlgorithmID = get(gui.AlgorithmMenu,'Value');
esfitdata.Opts.TargetID = get(gui.TargetMenu,'Value');
esfitdata.Opts.AutoScaleID = get(gui.AutoScaleMenu,'Value');
if esfitdata.nDataSets==1
  esfitdata.Opts.BaseLine = get(gui.BaseLineMenu,'Value');
else
  weight = zeros(1,sum(esfitdata.datasize));
  for i = 1:esfitdata.nDataSets
    esfitdata.Opts.BaseLine(i) = get(gui.multifitBaseline(i),'Value');
    esfitdata.Opts.weight(i) = str2double(get(gui.multifitWeight(i),'Value'));
    weight(esfitdata.idx{i}) = esfitdata.Opts.weight(i);
  end
  esfitdata.weight = weight;
end
esfitdata.Opts.useMask = get(gui.MaskCheckbox,'Value')==1;

% Run fitting
%-------------------------------------------------------------------------------
useGUI = true;
try
  result = runFitting(useGUI);
catch ME
  if isfield(esfitdata,'modelEvalError') && esfitdata.modelEvalError
    return
  elseif contains(ME.stack(1).name,'esfit')
    GUIErrorHandler(ME);
    return;
  else
    error(ME.message)
  end
end

% Save result to fit parameter set list
esfitdata.currFitSet = result;
esfitdata.currFitSet.Mask = esfitdata.Opts.useMask && ~all(esfitdata.Opts.Mask);

% Update GUI with fit results
%-------------------------------------------------------------------------------

% Remove current values and uncertainties from parameter table
Data = gui.ParameterTable.Data;
pi = 1;
for p = 1:size(Data,1)
  Data{p,7} = '-'; 
  if ~esfitdata.fixedParams(p) && ~isempty(esfitdata.best.pstd)
    Data{p,8} = sprintf('%0.6g',esfitdata.best.pfit(pi));
    Data{p,9} = sprintf('%0.6g',esfitdata.best.ci95(pi,1));
    Data{p,10} = sprintf('%0.6g',esfitdata.best.ci95(pi,2));
    pi = pi+1;
  else
    Data{p,9} = '-';
    Data{p,10} = '-';
  end
end
set(gui.ParameterTable,'Data',Data);

% Hide current sim plot in data axes
for i = 1:esfitdata.nDataSets
  idx = esfitdata.idx{i};
  expdata = esfitdata.data(idx);
  bestsim = esfitdata.best.fit(idx);
  set(gui.currsimdata(i),'YData',NaN(1,numel(expdata)));
  residualzero = get(gui.residualzero(i),'Value');
  set(gui.residualdata(i),'YData',(expdata(:)-bestsim(:))+residualzero);

  if esfitdata.nDataSets>1
    gui.multifitScale(i).Value = esfitdata.best.scale(i);
  end
end
drawnow

% Reactivate UI components
set(gui.SaveButton,'Enable','on');

if isfield(esfitdata,'FitSets') && numel(esfitdata.FitSets)>0
  set(gui.selectStartPointButtonSelected,'Enable','on');
  set(gui.exportSetButton,'Enable','on');
  set(gui.saveSetButton,'Enable','on');
  set(gui.deleteSetButton,'Enable','on');
  set(gui.clearSetButton,'Enable','on');
end

% Hide stop button, show start button
set(gui.StartButton,'Text','Start fitting');
set(gui.StartButton,'Tooltip','Start fitting');

% Re-enable other buttons
set(gui.EvaluateButton,'Enable','on');
set(gui.ResetButton,'Enable','on');
set(gui.AlgorithmSettingsButton,'Enable','on');

% Re-enable listboxes
set(gui.AlgorithmMenu,'Enable','on');
set(gui.TargetMenu,'Enable','on');
set(gui.AutoScaleMenu,'Enable','on');
set(gui.BaseLineMenu,'Enable','on');
if esfitdata.nDataSets>1
  set(gui.multifitBaseline,'Enable','on');
  set(gui.multifitWeight,'Enable','on');
end
set(gui.AutoScaleMenu,'Enable','on');

% Re-enable parameter table and its selection controls
set(gui.selectAllButton,'Enable','on');
set(gui.selectNoneButton,'Enable','on');
set(gui.selectInvButton,'Enable','on');
set(gui.selectStartPointButtonCenter,'Enable','on');
set(gui.selectStartPointButtonRandom,'Enable','on');
set(gui.selectStartPointButtonBest,'Enable','on');
set(gui.ParameterTable,'ColumnEditable',colEditable);
set(gui.ParameterTable,'CellEditCallback',@paramsTableCellEditCallback);

% Re-enable mask tools
for i = 1:numel(gui.dataaxes)
  gui.dataaxes(i).ButtonDownFcn = @(src,evt) axesButtonDownFcn(src);
end
set(gui.clearMaskButton,'Enable','on');
set(gui.MaskCheckbox,'Enable','on');

% Re-enable algorithm settings window
set(gui.AlgorithmMenuPopup,'Enable','on');
set(gui.setAlgorithmDefaults,'Enable','on');
nTabs = numel(gui.AlgorithmTabs.UserData);
for i = 1:nTabs
  hsetting = gui.AlgorithmTabs.UserData{i};
  for j = 1:numel(hsetting)
    set(hsetting(j),'Enable','on');
  end
end

% Re-enable multifit settings
if esfitdata.nDataSets>1
  % set(gui.TilesButton,'Enable','on');
  % set(gui.TabsButton,'Enable','on');
  set(gui.multifitWeight,'Enable','on');
  set(gui.multifitBaseline,'Enable','on');
end

end
%===============================================================================

%===============================================================================
function resetCallback(~,~)
% Reset best fit
global esfitdata gui

% Remove messages from log
scroll(gui.LogBox,'top')
set(gui.LogBox,'Value',{''})

% Remove simulations from plot
for i = 1:esfitdata.nDataSets

  gui.bestsimdata(i).YData = NaN(size(gui.bestsimdata(i).YData));
  gui.currsimdata(i).YData = NaN(size(gui.currsimdata(i).YData));
  gui.baselinedata(i).YData = NaN(size(gui.baselinedata(i).YData));
  gui.residualdata(i).YData = NaN(size(gui.residualdata(i).YData));

  if esfitdata.nDataSets>1
    gui.multifitScale(i).Value = 0;
  end

end

esfitdata.best = [];
esfitdata.curr = [];

% Remove mask patches
clearMaskCallback()
gui.MaskCheckbox.Value = 0;

% Readjust vertical range
updateaxislimits()

% Reset rmsdhistory plot
esfitdata.rmsdhistory = [];
updatermsdplot;
iterationprint('');

% Reset besthistory 
esfitdata.besthistory.rmsd = [];
esfitdata.besthistory.par = [];

% Remove displayed best fit and uncertainties
Data = gui.ParameterTable.Data;
for p = 1:size(Data,1)
  Data{p,2} = true;
  Data{p,7} = '-';
  Data{p,8} = '-';
  Data{p,9} = '-';
  Data{p,10} = '-';
end
set(gui.ParameterTable,'Data',Data);

% Reset algorithm defaults
updateAlgorithmDefaults();

% Reset plot
gui.BaselineCheckbox.Value = 0;
showbaseline(gui.BaselineCheckbox);
gui.ResidualCheckbox.Value = 0;
showresiduals(gui.ResidualCheckbox);

end
%===============================================================================

%===============================================================================
% Function to update GUI at each iteration
%-------------------------------------------------------------------------------
function userstop = iterupdateGUI(info)
global esfitdata gui
userstop = esfitdata.UserCommand~=0;
windowClosing = esfitdata.UserCommand==99;
if windowClosing, return; end

currpar = esfitdata.curr.par;
bestpar = esfitdata.best.par;

% Update plot(s)
for i = 1:esfitdata.nDataSets

  idx = esfitdata.idx{i};

  % Get relevant quantities
  x = esfitdata.Opts.x(idx);
  expdata = esfitdata.data(idx);
  bestsim = real(esfitdata.best.fit(idx));
  currsim = real(esfitdata.curr.sim(idx));
  currbaseline = real(esfitdata.curr.baseline(idx));

  % Update plotted data
  set(gui.expdata(i),'XData',x,'YData',expdata);
  set(gui.bestsimdata(i),'XData',x,'YData',bestsim);
  set(gui.currsimdata(i),'XData',x,'YData',currsim);
  set(gui.baselinedata(i),'XData',x,'YData',currbaseline);
  residualzero = get(gui.residualzero(i),'Value');
  residualdata = (expdata(:)-currsim(:))+residualzero;
  set(gui.residualdata(i),'XData',x,'YData',residualdata);

  if esfitdata.nDataSets>1
    gui.multifitScale(i).Value = esfitdata.curr.scale(i);
  end

end

% Readjust vertical range
updateaxislimits()

% Update column with current parameter values
data = get(gui.ParameterTable,'data');
nParams = size(data,1);
for p = 1:nParams
  oldvaluestring = striphtml(data{p,7});
  newvaluestring = sprintf('%0.6f',currpar(p));
  % Find first character at which the new value differs from the old one
  idx = 1;
  while idx<=min(length(oldvaluestring),length(newvaluestring))
    if oldvaluestring(idx)~=newvaluestring(idx), break; end
    idx = idx + 1;
  end
  if verLessThan('Matlab','9.13') % 9.13 = 2022b
    str = newvaluestring;
  else
    active = data{p,2};
    if active
      str = ['<html><font color="#000000">' newvaluestring(1:idx-1) '</font><font color="#888888">' newvaluestring(idx:end) '</font></html>'];
    else
      str = ['<html><font color="#888888">' newvaluestring '</font></html>'];
    end
    % Indicate parameters have hit limit
    if currpar(p)==esfitdata.pvec_lb(p) ||  currpar(p)==esfitdata.pvec_ub(p)
      str = ['<html><font color="#ff0000">' newvaluestring '</font></html>'];
    end
  end
  data{p,7} = str;
end

% Update column with best values if current parameter set is new best
if info.newbest

  str = sprintf('Current best RMSD: %g\n',esfitdata.best.rmsd);
  set(gui.RmsText,'Text',str,'FontColor',[0 0.6 0]);

  for p = 1:nParams
    oldvaluestring = striphtml(data{p,8});
    newvaluestring = sprintf('%0.6g',bestpar(p));
    % Find first character at which the new value differs from the old one
    idx = 1;
    while idx<=min(length(oldvaluestring),length(newvaluestring))
      if oldvaluestring(idx)~=newvaluestring(idx), break; end
      idx = idx + 1;
    end
    if verLessThan('Matlab','9.13') % 9.13 = 2022b
      str = newvaluestring;
    else
      active = data{p,2};
      if active
        if bestpar(p)==esfitdata.pvec_lb(p) ||  bestpar(p)==esfitdata.pvec_ub(p)
          str = ['<html><font color="#ff0000">' newvaluestring '</font></html>'];
        else
          str = ['<html><font color="#009900">' newvaluestring(1:idx-1) '</font><font color="#000000">' newvaluestring(idx:end) '</font></html>'];
        end
      else
        str = ['<html><font color="#888888">' newvaluestring '</font></html>'];
      end
    end
    data{p,8} = str;
  end

  if gui.autosave.Value

    % Assemble result structure
    fitresult.algorithm = esfitdata.AlgorithmAbbrev{esfitdata.Opts.AlgorithmID};
    fitresult.target = esfitdata.TargetAbbrev{esfitdata.Opts.TargetID};
    if esfitdata.nDataSets>1
      for k = 1:esfitdata.nDataSets
        idx = esfitdata.idx{k};
        fitresult.fit{k} = esfitdata.best.fit(idx);
        fitresult.fitraw{k} = esfitdata.best.fit(idx)/esfitdata.best.scale(k);
        fitresult.mask{k} = esfitdata.Opts.Mask(idx);
        fitresult.residuals{k} = esfitdata.best.residuals(idx);
        fitresult.baseline{k} = esfitdata.best.baseline(idx);
      end
    else
      fitresult.fit = esfitdata.best.fit;
      fitresult.fitraw = esfitdata.best.fit/esfitdata.best.scale;
      fitresult.mask = esfitdata.Opts.Mask;
      fitresult.residuals = esfitdata.best.residuals;
      fitresult.baseline = esfitdata.best.baseline;
    end

    fitresult.baselinetype = esfitdata.best.baselinetype;
    fitresult.scale = esfitdata.best.scale;

    fitresult.bestfithistory.rmsd = esfitdata.besthistory.rmsd;
    fitresult.bestfithistory.pfit = esfitdata.besthistory.par;
    if esfitdata.structureInputs
      fitresult.bestfithistory.pfit2structs = esfitdata.p2args;
    end

    fitresult.pnames = {esfitdata.pinfo.Name}.';
    fitresult.pnames = fitresult.pnames(~esfitdata.fixedParams);
    fitresult.p_start = esfitdata.p_start;
    fitresult.p_fixed = esfitdata.fixedParams;
    fitresult.pfit = esfitdata.best.par(~esfitdata.fixedParams);
    fitresult.pfit_full = esfitdata.best.par;

    if esfitdata.structureInputs
      argsfit = esfitdata.p2args(fitresult.pfit_full);
    else
      argsfit = [];
    end
    fitresult.argsfit = argsfit;

    % Autosave result structure
    path = pwd;
    file = 'esfit_autosave.mat';
    save(fullfile(path,file),'fitresult')
  end
end

gui.ParameterTable.Data = data;

updatermsdplot;

drawnow limitrate

end
%===============================================================================

%===============================================================================
function updatermsdplot(~,~)
global esfitdata gui
% Update rmsd plot
if isfield(esfitdata,'best') && ~isempty(esfitdata.best) && ~isempty(esfitdata.best.rmsd)
  str = sprintf('Current best RMSD: %g\n',esfitdata.best.rmsd);
else
  str = sprintf('Current best RMSD: -\n');
end
set(gui.RmsText,'Text',str);

if ~isempty(gui.rmsdline)
  n = min(100,numel(esfitdata.rmsdhistory));
  set(gui.rmsdline,'XData',1:n,'YData',esfitdata.rmsdhistory(end-n+1:end));
  ax = gui.rmsdline.Parent;
  axis(ax,'tight');
  if gui.RmsLogPlot.Value
    set(ax,'yscale','log')
  else
    set(ax,'yscale','linear')
  end
end
end
%===============================================================================

%===============================================================================
function iterationprint(str)
global gui
if isempty(gui)
  disp(strtrim(str));
else
  set(gui.logLine,'Text',strtrim(str));
end
end
%===============================================================================

%===============================================================================
function setListCallback(~,~)
displayFitSet
end
%===============================================================================

%===============================================================================
function displayFitSet
global esfitdata gui
if isempty(gui.FitSetTable.Selection), return; end
fitsets = gui.FitSetTable.Data;
if ~isempty(fitsets)
  k = find([esfitdata.FitSets.ID]==fitsets{gui.FitSetTable.Selection,1});
  if k>0
    fitset = esfitdata.FitSets(k);

    % Set column with best-fit parameter values
    data = get(gui.ParameterTable,'data');

    pi = 1;
    for p = 1:size(data,1)
      data{p,8} = sprintf('%0.6g',fitset.pfit_full(p));
      if ~fitset.p_fixed(p) && ~isempty(fitset.pstd)
        data{p,9} = sprintf('%0.6g',fitset.ci95(pi,1));
        data{p,10} = sprintf('%0.6g',fitset.ci95(pi,2));
        pi = pi+1;
      else
        data{p,9} = '-';
        data{p,10} = '-';
      end
    end
    set(gui.ParameterTable,'Data',data);

    if esfitdata.nDataSets==1
      set(gui.bestsimdata,'YData',fitset.fit);
      set(gui.baselinedata,'YData',fitset.baseline);

      expdata = get(gui.expdata,'YData');
      residualzero = get(gui.residualzero,'Value');
      residualdata = (expdata(:)-fitset.fit(:))+residualzero;
      set(gui.residualdata,'YData',residualdata);
    else
      for i = 1:esfitdata.nDataSets
        set(gui.bestsimdata(i),'YData',fitset.fit{i});
        set(gui.baselinedata(i),'YData',fitset.baseline{i});
        gui.multifitScale(i).Value = fitset.scale(i);

        expdata = get(gui.expdata(i),'YData');
        residualzero = get(gui.residualzero(i),'Value');
        residualdata = (expdata(:)-fitset.fit{i}(:))+residualzero;
        set(gui.residualdata(i),'YData',residualdata);
      end
    end
    drawnow limitrate
  end
end

end
%===============================================================================

%===============================================================================
function saveFitsetCallback(~,~)
global esfitdata
if ~isempty(esfitdata.currFitSet)
  esfitdata.lastSetID = esfitdata.lastSetID+1;
  esfitdata.currFitSet.ID = esfitdata.lastSetID;
  if ~isfield(esfitdata,'FitSets') || isempty(esfitdata.FitSets)
    esfitdata.FitSets(1) = esfitdata.currFitSet;
  else
    esfitdata.FitSets(end+1) = esfitdata.currFitSet;
  end
  refreshFitsetList(-1);
end
end
%===============================================================================

%===============================================================================
function refreshFitsetList(idx)
global esfitdata gui
nSets = numel(esfitdata.FitSets);
data = cell(nSets,4);
for k = 1:nSets
  data{k,1} = esfitdata.FitSets(k).ID;
  data{k,2} = esfitdata.FitSets(k).rmsd;
  data{k,3} = esfitdata.FitSets(k).algorithm;
  data{k,4} = esfitdata.FitSets(k).target;
  if esfitdata.nDataSets==1
    data{k,5} = esfitdata.FitSets(k).scale;
    data{k,6} = esfitdata.FitSets(k).baselinetype{1};
  else
    str = sprintf('%5.4g ',esfitdata.FitSets(k).scale);
    data{k,5} = strcat('[',str(1:end-1),']');
    str = sprintf('%s, ',esfitdata.FitSets(k).baselinetype{:});
    data{k,6} = strcat('[',str(1:end-2),']');
  end
  if esfitdata.FitSets(k).Mask
    maskstr = ' (mask)';
  else
    maskstr = '';
  end
  data{k,7} = maskstr;
end
set(gui.FitSetTable,'Data',data);
if idx>0, set(gui.FitSetTable,'Selection',idx); end
if idx==-1, set(gui.FitSetTable,'Selection',nSets); end % R2021b onwards

if nSets>0, state = 'on'; else, state = 'off'; end
set(gui.selectStartPointButtonSelected,'Enable',state);
set(gui.exportSetButton,'Enable',state);
set(gui.saveSetButton,'Enable',state);
set(gui.deleteSetButton,'Enable',state);
set(gui.clearSetButton,'Enable',state);
displayFitSet;
end
%===============================================================================

%===============================================================================
function exportSetButtonCallback(~,~)
global esfitdata gui
ID = gui.FitSetTable.Data{gui.FitSetTable.Selection,1};
idx = [esfitdata.FitSets.ID]==ID;
fitresult = esfitdata.FitSets(idx);
fitresult = rmfield(fitresult,'Mask');
varname = sprintf('fit%d',ID);
assignin('base',varname,fitresult);
fprintf('Fit parameter set %d assigned to variable ''%s''.\n',ID,varname);
evalin('base',varname);
end
%===============================================================================

%===============================================================================
function saveSetButtonCallback(~,~)
global esfitdata gui
ID = gui.FitSetTable.Data{gui.FitSetTable.Selection,1};
idx = [esfitdata.FitSets.ID]==ID;
fitresult = esfitdata.FitSets(idx);
fitresult = rmfield(fitresult,'Mask');
varname = sprintf('fit%d',ID);
[file,path] = uiputfile(strcat(varname,'.mat'));
if file~=0
  save(fullfile(path,file),'fitresult')
  figure(gui.Fig)
  updateLogBox(sprintf('Fit parameter set %d saved as %s in %s.\n',ID,file,path));
end
end
%===============================================================================

%===============================================================================
function deleteSetButtonCallback(~,~)
global esfitdata gui
idx = gui.FitSetTable.Selection;
nSets = numel(esfitdata.FitSets);
if nSets>0
    ID = gui.FitSetTable.Data{gui.FitSetTable.Selection,1};
    for k = nSets:-1:1
        if esfitdata.FitSets(k).ID==ID
            esfitdata.FitSets(k) = [];
        end
    end
    if idx>length(esfitdata.FitSets), idx = length(esfitdata.FitSets); end
    if idx==0, idx = 1; end
    gui.FitSetTable.Selection = idx;
    refreshFitsetList(0);
end

if isempty(gui.FitSetTable.Data)
    set(gui.selectStartPointButtonSelected,'Enable','off');
    set(gui.exportSetButton,'Enable','off');
    set(gui.saveSetButton,'Enable','off');
    set(gui.deleteSetButton,'Enable','off');
    set(gui.clearSetButton,'Enable','off');
end
end
%===============================================================================

%===============================================================================
function deleteSetListKeyPressFcn(src,event)
if strcmp(event.Key,'delete')
    deleteSetButtonCallback(src,event);
end
end
%===============================================================================

%===============================================================================
function clearSetButtonCallback(~,~)
global esfitdata gui
nSets = numel(esfitdata.FitSets);
for k = nSets:-1:1
  esfitdata.FitSets(k) = [];
end
refreshFitsetList(0);

esfitdata.lastSetID = 0;

set(gui.selectStartPointButtonSelected,'Enable','off');
set(gui.exportSetButton,'Enable','off');
set(gui.saveSetButton,'Enable','off');
set(gui.deleteSetButton,'Enable','off');
set(gui.clearSetButton,'Enable','off');
end
%===============================================================================

%===============================================================================
function setStartPoint(sel)

global esfitdata gui

% Update fixed parameters
Data = gui.ParameterTable.Data;
for p = 1:size(Data,1)
  esfitdata.fixedParams(p) = Data{p,2}==0;
end

% Set starting point
%-------------------------------------------------------------------------------
fixedParams = esfitdata.fixedParams;
activeParams = ~fixedParams;
nParameters = numel(esfitdata.pvec_0);
lb = esfitdata.pvec_lb;
ub = esfitdata.pvec_ub;
p_start = esfitdata.p_start;

switch sel
  case 'center'
    pcenter = (lb+ub)/2;
    p_start(activeParams) = pcenter(activeParams);
  case 'random' % random
    prandom = lb + rand(nParameters,1).*(ub-lb);
    p_start(activeParams) = prandom(activeParams);
  case 'selected'
    data = gui.FitSetTable.Data;
    if ~isempty(data)
      ID = data{gui.FitSetTable.Selection,1};
      idx = find([esfitdata.FitSets.ID]==ID);
      if ~isempty(idx)
        p_start = esfitdata.FitSets(idx).pfit_full;
      else
        error('Could not locate selected parameter set.');
      end
    else
      error('No saved parameter set yet.');
    end
  case 'best'
    if isfield(esfitdata,'best') && ~isempty(esfitdata.best)
      p_start(activeParams) = esfitdata.best.pfit;
    end
end
esfitdata.p_start = p_start;

% Check if new start values fall within bound range, adapt bounds if not
updatebounds = false;
if strcmp(sel,'selected') || strcmp(sel,'best')
  newlb = p_start<lb;
  if any(newlb)
    updatebounds = true;
    db = (ub(newlb)-lb(newlb))/2;
    esfitdata.pvec_lb(newlb) = p_start(newlb)-db;
  end
  newub = p_start>ub;
  if any(newub)
    updatebounds = true;
    db = (ub(newub)-lb(newub))/2;
    esfitdata.pvec_ub(newub) = p_start(newub)+db;
  end
  if updatebounds
    updateLogBox('Selected parameter set outside range. Adapting range.')
  end
end

% Update parameter table
data = get(gui.ParameterTable,'data');
for p = 1:numel(p_start)
  data{p,4} = sprintf('%0.6g',p_start(p));
  if updatebounds
    data{p,5} = sprintf('%0.6g',esfitdata.pvec_lb(p));
    data{p,6} = sprintf('%0.6g',esfitdata.pvec_ub(p));
  end
end
gui.ParameterTable.Data = data;

end
%===============================================================================

%===============================================================================
function selectButtonCallback(type)
global gui
data = gui.ParameterTable.Data;
switch type
  case 'all'
    data(:,2) = {true};
  case 'none'
    data(:,2) = {false};
  case 'invert'
    for k=1:size(data,1)
      data{k,2} = ~data{k,2};
    end
end

allParamsFixed = all(~cell2mat(data(:,2)));
activeParameters = sum(cell2mat(data(:,2)));
gui.nParamsField.Value = activeParameters;
gui.nParamsFieldPopup.Value = activeParameters;
if allParamsFixed
  set(gui.StartButton,'Enable','off');
else
  set(gui.StartButton,'Enable','on');
end

set(gui.ParameterTable,'Data',data);
end
%===============================================================================

%===============================================================================
function changeBaseLineSelection(src,~)
global gui
% Update baseline type selection for all datasets
if isequal(src,gui.BaseLineMenu)
  for i = 1:numel(gui.multifitBaseline)
    gui.multifitBaseline(i).Value = gui.BaseLineMenu.Value;
  end
else
  [uniquevalues,~,ic] = unique([gui.multifitBaseline.Value]);
  if numel(uniquevalues)==1
    gui.BaseLineMenu.Value = uniquevalues;
  else
    counts = accumarray(ic,1);
    [~,order] = sort(counts,'descend');
    gui.BaseLineMenu.Value = uniquevalues(order(1));
  end
end
end
%===============================================================================

%===============================================================================
function changemultifitWeight(src,evt)
try
  newValue = eval(src.Value);
  newValue = num2str(newValue);
catch
  newValue = evt.PreviousValue;
end
src.Value = newValue;
end
%===============================================================================

%===============================================================================
function showbaseline(src,~)
global gui
h = gui.baselinedata;
if src.Value
  set(h,'Visible','on')
else
  set(h,'Visible','off')
end
end
%===============================================================================

%===============================================================================
function showresiduals(src,~)

global gui

for i = 1:numel(gui.residualdata)
  if src.Value
    set(gui.residualdata(i),'Visible','on')
    set(gui.residualzero(i),'Visible','on')
  else
    set(gui.residualdata(i),'Visible','off')
    set(gui.residualzero(i),'Visible','off')
  end
end
updateaxislimits()

end
%===============================================================================

%===============================================================================
function updateaxislimits(~)

global esfitdata gui

for i = 1:esfitdata.nDataSets

  idx = esfitdata.idx{i};

  dataplots = gui.dataaxes(i).Children;

  % Get YData for all visible line plots in axes
  plottedData = [];
  for j = 1:numel(dataplots)
    if isa(dataplots(j),'matlab.graphics.primitive.Line') && ...
        dataplots(j).Visible && ~any(isnan(dataplots(j).YData))
      plottedData = [plottedData; dataplots(j).YData(:)];
      if strcmp(dataplots(j).DisplayName,'expdata')
        plottedExpData = dataplots(j).YData(:);
      end
    end
  end
  % Adjust limits based on experimental data and presence of residuals
  maxy = max(plottedExpData);
  miny = min(plottedExpData);
  if gui.ResidualCheckbox.Value
    miny = miny-2*abs((maxy-miny)/3);
  end
  % Check that all plotted data is within axis limits
  maxyall = max(plottedData);
  minyall = min(plottedData);
  if miny>minyall
    miny = minyall;
  end
  if maxy<maxyall
    maxy = maxyall;
  end
  YLimits = [miny maxy] + [-1 1]*esfitdata.Opts.PlotStretchFactor*(maxy-miny);

  set(gui.dataaxes(i),'YLim',YLimits);
  x = esfitdata.Opts.x(idx);
  minx = min(x);
  maxx = max(x);
  set(gui.dataaxes(i),'XLim',[minx maxx]);

  % Update tile number position
  if isfield(gui,'TilesButton') && gui.TilesButton.Value
    gui.tileID(i).Position(1:2) = [minx+0.05*(maxx-minx) YLimits(2)-0.05*diff(YLimits)];
  end

  % Readjust mask patches
  tmp = gui.dataaxes(i).UserData;
  if isfield(tmp,'maskPatch') && ~isempty(tmp.maskPatch)
    for mp = 1:numel(tmp.maskPatch)
      tmp.maskPatch(mp).YData = YLimits([1 1 2 2]).';
    end
  end

end

end
%===============================================================================


%===============================================================================
function paramsTableCellEditCallback(~,callbackData)
global esfitdata gui

% Get handle of table and row/column index of edited table cell
hTable = callbackData.Source;
ridx = callbackData.Indices(1);
cidx = callbackData.Indices(2);

if cidx==2
  allParamsFixed = all(~cell2mat(hTable.Data(:,2)));
  activeParameters = sum(cell2mat(hTable.Data(:,2)));
  gui.nParamsField.Value = activeParameters;
  gui.nParamsFieldPopup.Value = activeParameters;
  if allParamsFixed
    set(gui.StartButton,'Enable','off');
  else
    set(gui.StartButton,'Enable','on');
  end
end

% Return unless it's a cell that contains start value or lower or upper bound
startColumn = 4; % start value column
lbColumn = 5; % lower-bound column
ubColumn = 6; % upper-bound column
startedit = cidx==startColumn;
lbedit = cidx==lbColumn;
ubedit = cidx==ubColumn;
if ~startedit && ~lbedit && ~ubedit, return; end

% Convert user-entered string to number
newval = str2double(callbackData.EditData);

% Revert if conversion didn't yield a scalar
if numel(newval)~=1 || isnan(newval) || ~isreal(newval)
  updateLogBox(sprintf('Input ''%s'' is not a number. Reverting edit.',callbackData.EditData));
  hTable.Data{ridx,cidx} = callbackData.PreviousData;
  return
end

% Get start value, lower and upper bounds of interval from table
start = str2double(hTable.Data{ridx,startColumn});
lower = str2double(hTable.Data{ridx,lbColumn});
upper = str2double(hTable.Data{ridx,ubColumn});

% Set new lower/upper bound
if startedit
  start = newval;
  if start<lower || start>upper
    updateLogBox('Start value outside range. Adapting range.');
    if start<lower
      lower = start-0.5*(upper-lower);
      hTable.Data{ridx,lbColumn} = sprintf('%0.6g',lower);
    elseif start>upper
      upper = start+0.5*(upper-lower);
      hTable.Data{ridx,ubColumn} = sprintf('%0.6g',upper);
    end
  end
elseif lbedit
  lower = newval;
elseif ubedit
  upper = newval;
end

% Revert if lower bound would be above upper bound
if lower>upper
  updateLogBox('Lower bound is above upper bound. Reverting edit.');
  hTable.Data{ridx,cidx} = callbackData.PreviousData;
  return
end

% Adapt start value if it falls outside new range
updatestartvalue = false;
if lower>start
  start = lower;
  updatestartvalue = true;
end
if upper<start
  start = upper;
  updatestartvalue = true;
end
if updatestartvalue
  updateLogBox('Start value outside new range. Adapting start value.');
  hTable.Data{ridx,startColumn} = sprintf('%0.6g',start);
end

% Update start value, lower and upper bounds
esfitdata.p_start(ridx) = start;
esfitdata.pvec_lb(ridx) = lower;
esfitdata.pvec_ub(ridx) = upper;

end
%===============================================================================

%===============================================================================
function axesButtonDownFcn(src)
global esfitdata gui

dataaxes = src;
i = str2double(dataaxes.Tag);
idx = esfitdata.idx{i};
mask_i = esfitdata.Opts.Mask(idx);

% Get mouse-click point on axes
cp = dataaxes.CurrentPoint;
x = esfitdata.Opts.x(idx);
[~,ind] = min(abs(x-cp(1,1)));
x1 = x(ind);

% Create temporary patch updating with user mouse motion
maskColor = [1 1 1]*0.95;
tmpmask = patch(dataaxes,x1*ones(1,4),dataaxes.YLim([1 1 2 2]),maskColor,'Tag','maskPatch','EdgeColor','none');

% Move new patch to the back
c = dataaxes.Children([2:end 1]);
dataaxes.Children = c;

% Continuously update patch based on mouse position until next user click
for i = 1:numel(gui.dataaxes)
  zoom(gui.dataaxes(i),'off')
end
set(gui.Fig,'WindowButtonMotionFcn',@(hObject,eventdata) drawmaskedregion(dataaxes,tmpmask));
set(gui.Fig,'WindowButtonDownFcn',@(src,event)uiresume(src));
dataaxes.ButtonDownFcn = [];
uiwait(gui.Fig);
dataaxes.ButtonDownFcn = @(src,evt) axesButtonDownFcn(src);
set(gui.Fig,'WindowButtonMotionFcn',[])
set(gui.Fig,'WindowButtonDownFcn',[])

% Update masked regions
cp = dataaxes.CurrentPoint;
[~,ind] = min(abs(x-cp(1,1)));
x2 = x(ind);
maskrange = sort([x1 x2]);
mask_i(x>maskrange(1) & x<maskrange(2)) = 0;
esfitdata.Opts.Mask(idx) = mask_i;
delete(tmpmask);
showmaskedregions();
end
%===============================================================================

%===============================================================================
function drawmaskedregion(dataaxes,tmpmask)
cp = get(dataaxes,'CurrentPoint');
xdata = tmpmask.XData;
xdata(2:3) = cp(1,1);
set(tmpmask,'XData',xdata);
end
%===============================================================================

%===============================================================================
function showmaskedregions()

global esfitdata gui

for i = 1:esfitdata.nDataSets

  idx = esfitdata.idx{i};
  x = esfitdata.Opts.x(idx).';

  tmp = gui.dataaxes(i).UserData;

  % Delete existing mask patches
  if isfield(tmp,'maskPatch') && ~isempty(tmp.maskPatch)
    delete(tmp.maskPatch);
    tmp = rmfield(tmp,'maskPatch');
  end

  % Show masked-out regions
  maskColor = [1 1 1]*0.95;
  edges = find(diff([1; esfitdata.Opts.Mask(idx); 1]));
  excludedRegions = reshape(edges,2,[]).';
  excludedRegions(:,1) = excludedRegions(:,1)-1;
  upperlimit = numel(x);
  excludedRegions(excludedRegions>upperlimit) = upperlimit;
  excludedRegions = x(excludedRegions);

  % Add a patch for each masked region
  nMaskPatches = size(excludedRegions,1);
  for r = 1:nMaskPatches
    x_patch = excludedRegions(r,[1 2 2 1]);
    y_patch = gui.dataaxes(i).YLim([1 1 2 2]);
    tmp.maskPatch(r) = patch(gui.dataaxes(i),x_patch,y_patch,maskColor,'EdgeColor','none');
  end

  % Reorder so that mask patches are in the back
  c = gui.dataaxes(i).Children([nMaskPatches+1:end, 1:nMaskPatches]);
  gui.dataaxes(i).Children = c;

  gui.dataaxes(i).UserData = tmp;

end

drawnow limitrate

end
%===============================================================================

%===============================================================================
function clearMaskCallback(~,~)
global esfitdata

esfitdata.Opts.Mask = true(size(esfitdata.Opts.Mask));
showmaskedregions();
esfitdata.best = [];
esfitdata.rmsdhistory = [];
esfitdata.besthistory.rmsd = [];
esfitdata.besthistory.par = [];

% Readjust vertical range
updateaxislimits()
drawnow

end
%===============================================================================

%===============================================================================
function updateLogBox(msg)
global gui

txt = get(gui.LogBox,'Value');
if numel(txt)==1 && isempty(txt{1})
  txt = {};
end
if ~iscell(msg)
  msg = cellstr(msg);
end
% Highlight errors
iserror = false;
if any(contains(msg,'Simulation function error','IgnoreCase',true))
  iserror = true;
end
for i = 1:numel(msg)
  msg{i} = strrep(msg{i},'\n','');
  msgs = strsplit(msg{i},newline);
  for j = 1:numel(msgs)
%     if iserror
%       msgs{j} = ['<html><font color="#ff0000">' msgs{j} '</font></html>'];
%     end
    txt{end+1} = msgs{j};
  end
end
if iserror
  txt{end+1} = '';
end
set(gui.LogBox,'Value',txt)
drawnow limitrate
scroll(gui.LogBox,'bottom')

end
%===============================================================================

%===============================================================================
function copyLog(~,~)
% Copy log to clipboard
global gui
txt = get(gui.LogBox,'Value');
str = [];
for i = 1:numel(txt)
  row = sprintf('%s\t', txt{i});
  row(end) = newline;
  str = [str row];
end
clipboard('copy',str)
end
%===============================================================================

%===============================================================================
function GUIErrorHandler(ME)
global esfitdata gui

% Reactivate UI components
set(gui.SaveButton,'Enable','off');

if isfield(esfitdata,'FitSets') && numel(esfitdata.FitSets)>0
  set(gui.selectStartPointButtonSelected,'Enable','on');
  set(gui.exportSetButton,'Enable','on');
  set(gui.saveSetButton,'Enable','on');
  set(gui.deleteSetButton,'Enable','on');
  set(gui.clearSetButton,'Enable','on');
end

% Hide stop button, show start button
set(gui.StartButton,'Text','Start fitting');
set(gui.StartButton,'Tooltip','Start fitting');

% Re-enable other buttons
set(gui.EvaluateButton,'Enable','on');
set(gui.ResetButton,'Enable','on');

% Re-enable listboxes
set(gui.AlgorithmMenu,'Enable','on');
set(gui.TargetMenu,'Enable','on');
set(gui.AutoScaleMenu,'Enable','on');
set(gui.BaseLineMenu,'Enable','on');
if esfitdata.nDataSets>1
  set(gui.multifitBaseline,'Enable','on');
  set(gui.multifitWeight,'Enable','on');
end

% Re-enable parameter table and its selection controls
set(gui.selectAllButton,'Enable','on');
set(gui.selectNoneButton,'Enable','on');
set(gui.selectInvButton,'Enable','on');
set(gui.selectStartPointButtonCenter,'Enable','on');
set(gui.selectStartPointButtonRandom,'Enable','on');
set(gui.selectStartPointButtonBest,'Enable','on');
colEditable = get(gui.ParameterTable,'UserData');
set(gui.ParameterTable,'ColumnEditable',colEditable);
set(gui.ParameterTable,'CellEditCallback',@paramsTableCellEditCallback);

% Re-enable mask tools
set(gui.clearMaskButton,'Enable','on');
set(gui.MaskCheckbox,'Enable','on');

if contains(ME.stack(1).name,'esfit')
  updateLogBox({ME.stack(1).name,' error:',ME.message})
else
  updateLogBox({'Simulation function error:',ME.message})
end

end
%===============================================================================

%===============================================================================
function openAlgorithmSettings(~,~)
% Callback for opening (closing) algorithm settings popup window

global esfitdata gui

if strcmp(gui.Popup.Visible,'off')
  AlgorithmID = esfitdata.Opts.AlgorithmID;
  if AlgorithmID<=numel(gui.AlgorithmTabs.Children)
    gui.AlgorithmTabs.SelectedTab = gui.AlgorithmTabs.Children(AlgorithmID);
  end
  set(gui.nParamsFieldPopup,'Value',sum(~esfitdata.fixedParams));
  set(gui.Popup,'Visible','on')
elseif strcmp(gui.Popup.Visible,'on')
  set(gui.Popup,'Visible','off');
end

end
%===============================================================================

%===============================================================================
function selectAlgorithm(src)
% Callback for selection of fitting algorithm

global esfitdata gui

% Update selection in main GUI window and algorithm settings popup
if nargin==1
  AlgorithmID = src.Value;
  gui.AlgorithmMenu.Value = AlgorithmID;
  if AlgorithmID<=numel(gui.AlgorithmTabs.Children)
    gui.AlgorithmMenuPopup.Value = AlgorithmID;
  end

  if AlgorithmID<=numel(gui.AlgorithmTabs.Children)
    gui.AlgorithmTabs.SelectedTab = gui.AlgorithmTabs.Children(AlgorithmID);
  end
else
  AlgorithmID = gui.AlgorithmMenu.Value;
end
esfitdata.Opts.AlgorithmID = AlgorithmID;

% Update the FitOpt structure
N = gui.nParamsFieldPopup.Value; %#ok<NASGU> 

if AlgorithmID<=numel(gui.AlgorithmTabs.Children)
  hsetting = gui.AlgorithmTabs.UserData{AlgorithmID};
  for i = 1:numel(hsetting)
    setting = hsetting(i).UserData{1};
    switch hsetting(i).UserData{2}
      case 'num'
        esfitdata.Opts.(setting) = str2double(hsetting(i).Value);
      case 'logical'
        esfitdata.Opts.(setting) = hsetting(i).Value;
      case 'eval'
        esfitdata.Opts.(setting) = eval(hsetting(i).Value);
    end
  end
end

end
%===============================================================================

%===============================================================================
function updateAlgorithmDefaults()
% Populate fit options settings popup with default options for all algorithms

global esfitdata gui
nTabs = numel(gui.AlgorithmTabs.UserData);

nParams = gui.nParamsFieldPopup.Value;

for i = 1:nTabs
  hsetting = gui.AlgorithmTabs.UserData{i};
  FitOpt = esfitdata.AlgorithmDefaults{i};
  for j = 1:numel(hsetting)
    setting = hsetting(j).UserData{1};
    switch hsetting(j).UserData{2}
      case 'num'
        hsetting(j).Value = num2str(FitOpt.(setting),'%g');
      case 'logical'
        hsetting(j).Value = FitOpt.(setting);
      case 'eval'
        switch setting
          case 'SimplexPars'
            if nParams>2
              hsetting(j).Value = '[1, 1+2/N, 0.75-1/(2*N), 1-1/N]';
            else
              hsetting(j).Value = '[1, 2, 0.5, 0.5]';
            end
          case 'EliteCount'
            hsetting(j).Value = num2str(max(2,ceil(0.1*FitOpt.PopulationSize)),'%i');
          case 'nParticles'
            hsetting(j).Value = '20 + N*10';
          case 'SwarmParams'
            hsetting(j).Value = sprintf('[%g, %g, %g, %g]',FitOpt.(setting));
        end
    end
  end

end

% Update current FitOpt structure
selectAlgorithm();

end
%===============================================================================

%===============================================================================
function changeAlgorithmSetting(src,evt)
% Callback for changed fit options setting

global esfitdata gui
N = gui.nParamsFieldPopup.Value; %#ok<NASGU> % needed to eval() expressions

setting = src.UserData{1};

% Check that entered expression evaluates without errors
newvalue = evt.Value;
oldvalue = evt.PreviousValue;
if strcmp(src.UserData{2},'eval') || strcmp(src.UserData{2},'num')
  try
    tmp = eval(newvalue); 
    if strcmp(src.UserData{2},'num')
      newvalue = tmp;
      src.Value = num2str(newvalue,'%g');
    else
      src.Value = newvalue;
    end
  catch ME
    % reset to old value in case of errors
    updateLogBox({strrep(ME.message,'Error:',sprintf('Error setting algorithm option %s:',setting))});
    src.Value = oldvalue;
  end
end

% Update FitOpt structure for selected algorithm
AlgorithmID = get(gui.AlgorithmMenu,'Value');

IDselected = find(gui.AlgorithmTabs.Children==gui.AlgorithmTabs.SelectedTab);

if IDselected==AlgorithmID
  switch src.UserData{2}
    case 'num'
      esfitdata.Opts.(setting) = str2double(src.Value);
    case 'logical'
      esfitdata.Opts.(setting) = src.Value;
    case 'eval'
      esfitdata.Opts.(setting) = eval(src.Value);
  end
end

end
%===============================================================================

%===============================================================================
function resizeGUI(type)
% Callback to resize GUI and enforce minimum acceptable size

global gui

switch type
  case 'main'
    newPosition = gui.Fig.Position;
    minwidth = 900;
    minheight = 650;
  case 'popup'
    newPosition = gui.Popup.Position;
    minwidth = 300;
    minheight = 300;
end

if newPosition(3)<minwidth
  newPosition(3) = minwidth;
end
if newPosition(4)<minheight
  oldheight = newPosition(4);
  diff = oldheight-minheight;
  newPosition(2) = newPosition(2)+diff;
  newPosition(4) = minheight;
end

switch type
  case 'main'
    gui.Fig.Position = newPosition;
  case 'popup'
    gui.Popup.Position = newPosition;
end

end
%===============================================================================

%===============================================================================
function closeGUI(~,~)
% Callback to close GUI

global esfitdata
esfitdata.UserCommand = 99;
drawnow;
hFig = findall(0,'Tag','esfitFigure');
if ~isempty(hFig)
  delete(hFig);
end
hPopup = findall(0,'Tag','algorithmpopup');
if ~isempty(hPopup)
  delete(hPopup);
end

end
%===============================================================================

%===============================================================================
function str = striphtml(str)
html = false;
for k = 1:numel(str)
  if ~html
    rmv(k) = false;
    if str(k)=='<', html = true; rmv(k) = true; end
  else
    rmv(k) = true;
    if str(k)=='>', html = false; end
  end
end
str(rmv) = [];
end
%===============================================================================

