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
%        .weights array of weights to use when combining residual vectors
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

function result = esfit_legacy(data,fcn,p0,varargin)

if nargin==0, help(mfilename); return; end

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

TargetNames{1} = 'data as is';
TargetNames{2} = 'integral';
TargetNames{3} = 'double integral';
TargetNames{4} = 'derivative';
TargetNames{5} = 'Fourier transform';
esfitdata.TargetNames = TargetNames;

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
  if EasySpinFunction
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

% Weights for global fitting
if ~isfield(Opt,'weights')
  Opt.weights = ones(1,esfitdata.nDataSets);
end
if numel(Opt.weights)~=esfitdata.nDataSets
  error('The number of elements in Opt.weights must be equal to the number of datasets.');
end

weights = zeros(1,sum(esfitdata.datasize));
for i = 1:numel(data)
  weights(esfitdata.idx{i}) = Opt.weights(i);
end
esfitdata.weights = weights;

% x axis for plotting
esfitdata.showxaxis = false;
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
  esfitdata.showxaxis = true;
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

% Setup GUI and return if in interactive mode
%-------------------------------------------------------------------------------
interactiveMode = nargout==0;
Opt.InfoPrintFunction = @(str) infoprint(str,interactiveMode);
esfitdata.Opts = Opt;
if interactiveMode
  if esfitdata.nDataSets>1
    error('The esfit GUI currently does not support global fitting.')
  end
  setupGUI(data);
  return
end

% Close GUI window if running in script mode
%-------------------------------------------------------------------------------
hFig = findobj('Tag','esfitFigure');
if ~isempty(hFig)
  close(hFig);
  esfitdata.UserCommand = 0;
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

end

%===============================================================================
%===============================================================================
%===============================================================================


%===============================================================================
% Run fitting algorithm
%===============================================================================
function result = runFitting(useGUI)

if nargin<1, useGUI = false; end

global esfitdata
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
  fitraw(idx) = fit(idx)/scale(i);
end
baseline = esfitdata.best.baseline;

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
calculateUncertainties = esfitdata.UserCommand==0 && nActiveParams>0;
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
    msg = 'Fitting stopped by user. Skipping uncertainty quantification.';
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
  result.fit = fit;
  result.fitraw = fitraw;
  result.mask = esfitdata.Opts.Mask;
  result.residuals = residuals0;
  result.baseline = baseline;
end

result.baseline = baseline;
result.scale = scale;
result.mask = esfitdata.Opts.Mask;

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

result.residuals = residuals0;
result.ssr = ssr0;
result.rmsd = rmsd0;
result.redchisquare = red_chisquare;

esfitdata.best = result;

end



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

simdata = out{esfitdata.OutArgument(1)}; % pick appropriate output argument
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

% Get x-axis if additional output argument is defined
if numel(esfitdata.OutArgument)==2
  xdata = out{esfitdata.OutArgument(2)};
  if iscell(xdata)
    x_vec = [];
    for i = 1:numel(xdata)
      x_vec = [x_vec; xdata{i}(:)];
    end
  else
    x_vec = xdata;
  end
  if numel(x_vec)~=numel(simdata_vec)
    error('\n  Simulation function output axis and data have unequal total number of points:\n    axis: %d\n    data: %d\n',...
           numel(x_vec),numel(simdata_vec));
  end
  esfitdata.Opts.x = x_vec;
end

% Rescale simulated data if scale should be ignored; include baseline if wanted
%-------------------------------------------------------------------------------
order = Opt.BaseLine;
baseline = zeros(size(simdata_vec));
simscale = ones(1,esfitdata.nDataSets);
for k = 1:esfitdata.nDataSets
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
rmsd = sqrt(mean(abs(residuals).^2.*esfitdata.weights(:)));
rmsd0 = sqrt(mean(abs(residuals0).^2));

esfitdata.curr.sim = simdata_vec;
esfitdata.curr.par = par;
esfitdata.curr.scale = simscale;
esfitdata.curr.baseline = baseline;

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
% Function to update GUI at each iteration
%-------------------------------------------------------------------------------
function userstop = iterupdateGUI(info)
global esfitdata
userstop = esfitdata.UserCommand~=0;
windowClosing = esfitdata.UserCommand==99;
if windowClosing, return; end

% Get relevant quantities
x = esfitdata.Opts.x(:);
expdata = esfitdata.data(:);
bestsim = real(esfitdata.best.fit(:));
currsim = real(esfitdata.curr.sim(:));
currpar = esfitdata.curr.par;
bestpar = esfitdata.best.par;

% Update plotted data
set(findobj('Tag','expdata'),'XData',x,'YData',expdata);
set(findobj('Tag','bestsimdata'),'XData',x,'YData',bestsim);
set(findobj('Tag','currsimdata'),'XData',x,'YData',currsim);

% Readjust vertical range
mask = esfitdata.Opts.Mask;
plottedData = [expdata(mask); bestsim; currsim];
maxy = max(plottedData);
miny = min(plottedData);
YLimits = [miny maxy] + [-1 1]*esfitdata.Opts.PlotStretchFactor*(maxy-miny);
hAx = findobj('Tag','dataaxes');
set(hAx,'YLim',YLimits);
set(hAx,'XLim',[min(x) max(x)]);

% Readjust mask patches
maskPatches = findobj('Tag','maskPatch');
for mp = 1:numel(maskPatches)
  maskPatches(mp).YData = YLimits([1 1 2 2]).';
end

% Update column with current parameter values
hParamTable = findobj('Tag','ParameterTable');
data = get(hParamTable,'data');
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
  data{p,7} = str;
end

% Update column with best values if current parameter set is new best
if info.newbest

  str = sprintf('Current best RMSD: %g\n',esfitdata.best.rmsd);
  hRmsText = findobj('Tag','RmsText');
  set(hRmsText,'String',str,'ForegroundColor',[0 0.6 0]);

  for p = 1:nParams
    oldvaluestring = striphtml(data{p,8});
    newvaluestring = sprintf('%0.6g',bestpar(p));
    % Find first character at which the new value differs from the old one
    idx = 1;
    while idx<=min(length(oldvaluestring),length(newvaluestring))
      if oldvaluestring(idx)~=newvaluestring(idx), break; end
      idx = idx + 1;
    end
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
    data{p,8} = str;
  end
end

hParamTable.Data = data;

updatermsdplot;

drawnow

end
%===============================================================================

%===============================================================================
function updatermsdplot(~,~)
global esfitdata
% Update rmsd plot
hRmsText = findobj('Tag','RmsText');
if isfield(esfitdata,'best') && ~isempty(esfitdata.best) && ~isempty(esfitdata.best.rmsd)
  str = sprintf('Current best RMSD: %g\n',esfitdata.best.rmsd);
else
  str = sprintf('Current best RMSD: -\n');
end
set(hRmsText,'String',str);

hRmsLogPlot = findobj('Tag','RmsLogPlot');

hrmsdline = findobj('Tag','rmsdline');
if ~isempty(hrmsdline)
  n = min(100,numel(esfitdata.rmsdhistory));
  set(hrmsdline,'XData',1:n,'YData',esfitdata.rmsdhistory(end-n+1:end));
  ax = hrmsdline.Parent;
  axis(ax,'tight');
  if hRmsLogPlot.Value
    set(ax,'yscale','log')
  else
    set(ax,'yscale','linear')
  end
end
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
function startButtonCallback(~,~)

global esfitdata
esfitdata.UserCommand = 0;

% Update GUI
%-------------------------------------------------------------------------------
% Hide Start button, show Stop button
set(findobj('Tag','StopButton'),'Visible','on');
set(findobj('Tag','StartButton'),'Visible','off');
set(findobj('Tag','SaveButton'),'Enable','off');

% Disable other buttons
set(findobj('Tag','EvaluateButton'),'Enable','off');
set(findobj('Tag','ResetButton'),'Enable','off');

% Disable listboxes
set(findobj('Tag','AlgorithmMenu'),'Enable','off');
set(findobj('Tag','TargetMenu'),'Enable','off');
set(findobj('Tag','BaseLineMenu'),'Enable','off');
set(findobj('Tag','AutoScaleMenu'),'Enable','off');

% Disable parameter table
set(findobj('Tag','selectAllButton'),'Enable','off');
set(findobj('Tag','selectNoneButton'),'Enable','off');
set(findobj('Tag','selectInvButton'),'Enable','off');
set(findobj('Tag','selectStartPointButtonCenter'),'Enable','off');
set(findobj('Tag','selectStartPointButtonRandom'),'Enable','off');
set(findobj('Tag','selectStartPointButtonSelected'),'Enable','off');
set(findobj('Tag','selectStartPointButtonBest'),'Enable','off');
colEditable = get(findobj('Tag','ParameterTable'),'UserData');
set(findobj('Tag','ParameterTable'),'ColumnEditable',false(size(colEditable)));
set(findobj('Tag','ParameterTable'),'CellEditCallback',[]);

% Remove displayed best fit and uncertainties
hTable = findobj('Tag','ParameterTable');
Data = hTable.Data;
for p = 1:size(Data,1)
  Data{p,7} = '-';
  Data{p,8} = '-';
  Data{p,9} = '-';
  Data{p,10} = '-';
end
set(hTable,'Data',Data);

% Get fixed parameters
for p = 1:esfitdata.nParameters
  esfitdata.fixedParams(p) = Data{p,2}==0;
end

% Disable fitset list controls
set(findobj('Tag','deleteSetButton'),'Enable','off');
set(findobj('Tag','exportSetButton'),'Enable','off');
set(findobj('Tag','sortIDSetButton'),'Enable','off');
set(findobj('Tag','sortRMSDSetButton'),'Enable','off');

% Disable mask tools
hAx = findobj('Tag','dataaxes');
hAx.ButtonDownFcn = [];
set(findobj('Tag','clearMaskButton'),'Enable','off');
set(findobj('Tag','MaskCheckbox'),'Enable','off');

% Pull settings from UI
%-------------------------------------------------------------------------------
% Determine selected method, target, autoscaling, start point
esfitdata.Opts.AlgorithmID = get(findobj('Tag','AlgorithmMenu'),'Value');
esfitdata.Opts.TargetID = get(findobj('Tag','TargetMenu'),'Value');
esfitdata.Opts.AutoScaleID = esfitdata.AutoScaleSettings{get(findobj('Tag','AutoScaleMenu'),'Value')};
esfitdata.Opts.BaseLine = esfitdata.BaseLineSettings{get(findobj('Tag','BaseLineMenu'),'Value')};
esfitdata.Opts.useMask = get(findobj('Tag','MaskCheckbox'),'Value')==1;

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

% Save result to fit set list
esfitdata.currFitSet = result;
esfitdata.currFitSet.Mask = esfitdata.Opts.useMask && ~all(esfitdata.Opts.Mask);


% Update GUI with fit results
%-------------------------------------------------------------------------------

% Remove current values and uncertainties from parameter table
hTable = findobj('Tag','ParameterTable');
Data = hTable.Data;
pi = 1;
for p = 1:size(Data,1)
  Data{p,7} = '-'; 
  if ~esfitdata.fixedParams(p) && ~isempty(esfitdata.best.pstd)
    Data{p,9} = sprintf('%0.6g',esfitdata.best.ci95(pi,1));
    Data{p,10} = sprintf('%0.6g',esfitdata.best.ci95(pi,2));
    pi = pi+1;
  else
    Data{p,9} = '-';
    Data{p,10} = '-';
  end
end
set(hTable,'Data',Data);

% Hide current sim plot in data axes
set(findobj('Tag','currsimdata'),'YData',NaN(1,numel(esfitdata.data)));
drawnow

% Reactivate UI components
set(findobj('Tag','SaveButton'),'Enable','on');

if isfield(esfitdata,'FitSets') && numel(esfitdata.FitSets)>0
  set(findobj('Tag','deleteSetButton'),'Enable','on');
  set(findobj('Tag','exportSetButton'),'Enable','on');
  set(findobj('Tag','sortIDSetButton'),'Enable','on');
  set(findobj('Tag','sortRMSDSetButton'),'Enable','on');
end

% Hide stop button, show start button
set(findobj('Tag','StopButton'),'Visible','off');
set(findobj('Tag','StartButton'),'Visible','on');

% Re-enable other buttons
set(findobj('Tag','EvaluateButton'),'Enable','on');
set(findobj('Tag','ResetButton'),'Enable','on');

% Re-enable listboxes
set(findobj('Tag','AlgorithmMenu'),'Enable','on');
set(findobj('Tag','TargetMenu'),'Enable','on');
set(findobj('Tag','BaseLineMenu'),'Enable','on');
set(findobj('Tag','AutoScaleMenu'),'Enable','on');

% Re-enable parameter table and its selection controls
set(findobj('Tag','selectAllButton'),'Enable','on');
set(findobj('Tag','selectNoneButton'),'Enable','on');
set(findobj('Tag','selectInvButton'),'Enable','on');
set(findobj('Tag','selectStartPointButtonCenter'),'Enable','on');
set(findobj('Tag','selectStartPointButtonRandom'),'Enable','on');
set(findobj('Tag','selectStartPointButtonSelected'),'Enable','on');
set(findobj('Tag','selectStartPointButtonBest'),'Enable','on');
set(findobj('Tag','ParameterTable'),'ColumnEditable',colEditable);
set(findobj('Tag','ParameterTable'),'CellEditCallback',@tableCellEditCallback);

% Re-enable mask tools
hAx = findobj('Tag','dataaxes');
hAx.ButtonDownFcn = @axesButtonDownFcn;
set(findobj('Tag','clearMaskButton'),'Enable','on');
set(findobj('Tag','MaskCheckbox'),'Enable','on');

end
%===============================================================================


%===============================================================================
function iterationprint(str)
hLogLine = findobj('Tag','logLine');
if isempty(hLogLine)
  disp(strtrim(str));
else
  set(hLogLine,'String',strtrim(str));
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


%===============================================================================
function evaluateCallback(~,~)
% Evaluate for selected parameters
global esfitdata

esfitdata.modelErrorHandler = @(ME) GUIErrorHandler(ME);

p_eval = esfitdata.p_start;
active = ~esfitdata.fixedParams;
p_eval = p_eval(active);
expdata = esfitdata.data(:);
esfitdata.Opts.AutoScaleID = esfitdata.AutoScaleSettings{get(findobj('Tag','AutoScaleMenu'),'Value')};
esfitdata.Opts.BaseLine = esfitdata.BaseLineSettings{get(findobj('Tag','BaseLineMenu'),'Value')};
esfitdata.Opts.useMask = get(findobj('Tag','MaskCheckbox'),'Value')==1;
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

% Get current spectrum
currsim = real(esfitdata.curr.sim(:));

% Update plotted data
x = esfitdata.Opts.x(:);
set(findobj('Tag','expdata'),'XData',x,'YData',expdata);
set(findobj('Tag','currsimdata'),'XData',x,'YData',currsim);

% Readjust vertical range
mask = esfitdata.Opts.Mask;
if isfield(esfitdata,'best') && isfield(esfitdata.best,'fit')
  bestsim = real(esfitdata.best.fit(:));
else
  bestsim = zeros(size(currsim));
end
plottedData = [expdata(mask); bestsim; currsim];
maxy = max(plottedData);
miny = min(plottedData);
YLimits = [miny maxy] + [-1 1]*esfitdata.Opts.PlotStretchFactor*(maxy-miny);
hAx = findobj('Tag','dataaxes');
set(hAx,'YLim',YLimits);
set(hAx,'XLim',[min(x) max(x)]);
drawnow

% Readjust mask patches
maskPatches = findobj('Tag','maskPatch');
for mp = 1:numel(maskPatches)
  maskPatches(mp).YData = YLimits([1 1 2 2]).';
end

% Update column with best values if current parameter set is new best
str = sprintf('Current RMSD: %g\n',rmsd);
hRmsText = findobj('Tag','RmsText');
set(hRmsText,'String',str,'ForegroundColor',[1 0 0]);

end
%===============================================================================

%===============================================================================
function resetCallback(~,~)
% Reset best fit
global esfitdata

% Remove messages from log
set(findobj('Tag','LogBox'),'ListboxTop',1)
set(findobj('Tag','LogBox'),'String',{''})

% Remove best fit simulation from plot
hBestSim = findobj('Tag','bestsimdata');
hBestSim.YData = NaN(size(hBestSim.YData));
esfitdata.best = [];

% Remove current fit simulation from plot
hCurrSim = findobj('Tag','currsimdata');
hCurrSim.YData = NaN(size(hBestSim.YData));
esfitdata.curr = [];

% Readjust vertical range
mask = esfitdata.Opts.Mask;
expdata = esfitdata.data(:);
maxy = max(expdata(mask));
miny = min(expdata(mask));
YLimits = [miny maxy] + [-1 1]*esfitdata.Opts.PlotStretchFactor*(maxy-miny);
set(findobj('Tag','dataaxes'),'YLim',YLimits);
drawnow

% Reset rmsdhistory plot
esfitdata.rmsdhistory = [];
updatermsdplot;
iterationprint('');

% Reset besthistory 
esfitdata.besthistory.rmsd = [];
esfitdata.besthistory.par = [];

% Remove displayed best fit and uncertainties
hTable = findobj('Tag','ParameterTable');
Data = hTable.Data;
for p = 1:size(Data,1)
  Data{p,7} = '-';
  Data{p,8} = '-';
  Data{p,9} = '-';
  Data{p,10} = '-';
end
set(hTable,'Data',Data);

end
%===============================================================================

%===============================================================================
function deleteSetButtonCallback(~,~)
global esfitdata
h = findobj('Tag','SetListBox');
idx = h.Value;
str = h.String;
nSets = numel(str);
if nSets>0
    ID = sscanf(str{idx},'%d');
    for k = numel(esfitdata.FitSets):-1:1
        if esfitdata.FitSets(k).ID==ID
            esfitdata.FitSets(k) = [];
        end
    end
    if idx>length(esfitdata.FitSets), idx = length(esfitdata.FitSets); end
    if idx==0, idx = 1; end
    h.Value = idx;
    refreshFitsetList(0);
end

str = h.String';
if isempty(str)
    set(findobj('Tag','deleteSetButton'),'Enable','off');
    set(findobj('Tag','exportSetButton'),'Enable','off');
    set(findobj('Tag','sortIDSetButton'),'Enable','off');
    set(findobj('Tag','sortRMSDSetButton'),'Enable','off');
end
end
%===============================================================================


%===============================================================================
function deleteSetListKeyPressFcn(src,event)
if strcmp(event.Key,'delete')
    deleteSetButtonCallback(src,gco,event);
    displayFitSet
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
global esfitdata
h = findobj('Tag','SetListBox');
idx = h.Value;
str = h.String;
if ~isempty(str)
  ID = sscanf(str{idx},'%d');
  k = find([esfitdata.FitSets.ID]==ID);
  if k>0
    fitset = esfitdata.FitSets(k);

    % Set column with best-fit parameter values
    hTable = findobj('Tag','ParameterTable');
    data = get(hTable,'data');

    pi = 1;
    for p = 1:size(data,1)
      data{p,8} = sprintf('%0.6g',fitset.pfit_full(p));
      if ~fitset.fixedParams(p) && ~isempty(fitset.pstd)
        data{p,9} = sprintf('%0.6g',fitset.ci95(pi,1));
        data{p,10} = sprintf('%0.6g',fitset.ci95(pi,2));
        pi = pi+1;
      else
        data{p,9} = '-';
        data{p,10} = '-';
      end
    end
    set(hTable,'Data',data);

    h = findobj('Tag','bestsimdata');
    set(h,'YData',fitset.fit);
    drawnow
  end
end

end
%===============================================================================


%===============================================================================
function exportSetButtonCallback(~,~)
global esfitdata
h = findobj('Tag','SetListBox');
v = h.Value;
s = h.String;
ID = sscanf(s{v},'%d');
idx = [esfitdata.FitSets.ID]==ID;
fitresult = esfitdata.FitSets(idx);
fitresult = rmfield(fitresult,'Mask');
varname = sprintf('fit%d',ID);
assignin('base',varname,fitresult);
fprintf('Fit set %d assigned to variable ''%s''.\n',ID,varname);
evalin('base',varname);
end
%===============================================================================


%===============================================================================
function selectAllButtonCallback(~,~)
h = findobj('Tag','ParameterTable');
d = h.Data;
d(:,2) = {true};
set(h,'Data',d);
end
%===============================================================================


%===============================================================================
function selectNoneButtonCallback(~,~)
h = findobj('Tag','ParameterTable');
d = h.Data;
d(:,2) = {false};
set(h,'Data',d);
end
%===============================================================================


%===============================================================================
function selectInvButtonCallback(~,~)
h = findobj('Tag','ParameterTable');
d = h.Data;
for k=1:size(d,1)
    d{k,2} = ~d{k,2};
end
set(h,'Data',d);
end
%===============================================================================


%===============================================================================
function sortIDSetButtonCallback(~,~)
global esfitdata
ID = [esfitdata.FitSets.ID];
[~,idx] = sort(ID);
esfitdata.FitSets = esfitdata.FitSets(idx);
refreshFitsetList(0);
end
%===============================================================================


%===============================================================================
function sortRMSDSetButtonCallback(~,~)
global esfitdata
rmsd = [esfitdata.FitSets.rmsd];
[~,idx] = sort(rmsd);
esfitdata.FitSets = esfitdata.FitSets(idx);
refreshFitsetList(0);
end
%===============================================================================


%===============================================================================
function refreshFitsetList(idx)
global esfitdata
h = findobj('Tag','SetListBox');
nSets = numel(esfitdata.FitSets);
s = cell(1,nSets);
for k = 1:nSets
  if esfitdata.FitSets(k).Mask
    maskstr = ' (mask)';
  else
    maskstr = '';
  end
  s{k} = sprintf('%d. rmsd %g%s',...
    esfitdata.FitSets(k).ID,esfitdata.FitSets(k).rmsd,maskstr);
end
set(h,'String',s);
if idx>0, set(h,'Value',idx); end
if idx==-1, set(h,'Value',numel(s)); end

if nSets>0, state = 'on'; else, state = 'off'; end
set(findobj('Tag','deleteSetButton'),'Enable',state);
set(findobj('Tag','exportSetButton'),'Enable',state);
set(findobj('Tag','sortIDSetButton'),'Enable',state);
set(findobj('Tag','sortRMSDSetButton'),'Enable',state);

displayFitSet;
end
%===============================================================================

%===============================================================================
function GUIErrorHandler(ME)
global esfitdata

% Reactivate UI components
set(findobj('Tag','SaveButton'),'Enable','off');

if isfield(esfitdata,'FitSets') && numel(esfitdata.FitSets)>0
  set(findobj('Tag','deleteSetButton'),'Enable','on');
  set(findobj('Tag','exportSetButton'),'Enable','on');
  set(findobj('Tag','sortIDSetButton'),'Enable','on');
  set(findobj('Tag','sortRMSDSetButton'),'Enable','on');
end

% Hide stop button, show start button
set(findobj('Tag','StopButton'),'Visible','off');
set(findobj('Tag','StartButton'),'Visible','on');

% Re-enable other buttons
set(findobj('Tag','EvaluateButton'),'Enable','on');
set(findobj('Tag','ResetButton'),'Enable','on');

% Re-enable listboxes
set(findobj('Tag','AlgorithmMenu'),'Enable','on');
set(findobj('Tag','TargetMenu'),'Enable','on');
set(findobj('Tag','BaseLineMenu'),'Enable','on');
set(findobj('Tag','AutoScaleMenu'),'Enable','on');

% Re-enable parameter table and its selection controls
set(findobj('Tag','selectAllButton'),'Enable','on');
set(findobj('Tag','selectNoneButton'),'Enable','on');
set(findobj('Tag','selectInvButton'),'Enable','on');
set(findobj('Tag','selectStartPointButtonCenter'),'Enable','on');
set(findobj('Tag','selectStartPointButtonRandom'),'Enable','on');
set(findobj('Tag','selectStartPointButtonSelected'),'Enable','on');
set(findobj('Tag','selectStartPointButtonBest'),'Enable','on');
colEditable = get(findobj('Tag','ParameterTable'),'UserData');
set(findobj('Tag','ParameterTable'),'ColumnEditable',colEditable);
set(findobj('Tag','ParameterTable'),'CellEditCallback',@tableCellEditCallback);

% Re-enable mask tools
set(findobj('Tag','clearMaskButton'),'Enable','on');
set(findobj('Tag','MaskCheckbox'),'Enable','on');

if contains(ME.stack(1).name,'esfit')
  updateLogBox({ME.stack(1).name,' error:',ME.message})
else
  updateLogBox({'Simulation function error:',ME.message})
end

end
%===============================================================================

%===============================================================================
function updateLogBox(msg)

txt = get(findobj('Tag','LogBox'),'String');
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
    if iserror
      msgs{j} = ['<html><font color="#ff0000">' msgs{j} '</font></html>'];
    end
    txt{end+1} = msgs{j};
  end
end
if iserror
  txt{end+1} = '';
end
nval = numel(txt);
set(findobj('Tag','LogBox'),'String',txt)
drawnow
if nval>6
  txt = get(findobj('Tag','LogBox'),'String');
  set(findobj('Tag','LogBox'),'ListBoxTop',numel(txt)-5);
end

end
%===============================================================================

%===============================================================================
function copyLog(~,~)
% Copy log to clipboard
txt = get(findobj('Tag','LogBox'),'String');
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
function clearMaskCallback(~,~)
global esfitdata
esfitdata.Opts.Mask = true(size(esfitdata.Opts.Mask));
showmaskedregions();
esfitdata.best = [];
esfitdata.rmsdhistory = [];
esfitdata.besthistory.rmsd = [];
esfitdata.besthistory.par = [];

% Readjust vertical range
mask = esfitdata.Opts.Mask;
expdata = esfitdata.data(:);
maxy = max(expdata(mask));
miny = min(expdata(mask));
YLimits = [miny maxy] + [-1 1]*esfitdata.Opts.PlotStretchFactor*(maxy-miny);
set(findobj('Tag','dataaxes'),'YLim',YLimits);
drawnow

end
%===============================================================================


%===============================================================================
function saveFitsetCallback(~,~)
global esfitdata
if ~isempty(esfitdata.currFitSet)
  esfitdata.lastSetID = esfitdata.lastSetID+1;
  esfitdata.currFitSet.ID = esfitdata.lastSetID;
  esfitdata.currFitSet.fixedParams = esfitdata.fixedParams;
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
function tableCellEditCallback(~,callbackData)
global esfitdata

% Get handle of table and row/column index of edited table cell
hTable = callbackData.Source;
ridx = callbackData.Indices(1);
cidx = callbackData.Indices(2);

if cidx==1
  allParamsFixed = all(~cell2mat(hTable.Data(:,1)));
  if allParamsFixed
    set(findobj('Tag','StartButton'),'Enable','off');
  else
    set(findobj('Tag','StartButton'),'Enable','on');
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
    updateLogBox('Start value outside range. Reverting edit.');
    hTable.Data{ridx,cidx} = callbackData.PreviousData;
    return
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
function setupGUI(data)

if iscell(data)
  data = data{1};
end

global esfitdata
Opt = esfitdata.Opts;

% Main figure
%-------------------------------------------------------------------------------
hFig = findobj('Tag','esfitFigure');
if isempty(hFig)
  hFig = figure('Tag','esfitFigure','WindowStyle','normal');
else
  figure(hFig);
  clf(hFig);
end
set(hFig,'Visible','off')

sz = [1330 800]; % figure size
screensize = get(0,'ScreenSize');
scalefact = min(0.9*(screensize(3:4)/sz));
if scalefact>1
  scalefact = 1;
end
sz = sz*scalefact;
xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the screen vertically
set(hFig,'position',[xpos, ypos, sz(1), sz(2)],'units','pixels');
set(hFig,'WindowStyle','normal','DockControls','off','MenuBar','none');
set(hFig,'Resize','off');
set(hFig,'Name','esfit - Least-Squares Fitting','NumberTitle','off');
set(hFig,'CloseRequestFcn',...
    'global esfitdata; esfitdata.UserCommand = 99; drawnow; delete(gcf);');
  
spacing = 30*scalefact;
hPtop = 180*scalefact;
wPright = 230*scalefact;

Axesw = sz(1)-2*spacing-wPright;
Axesh = sz(2)-2.5*spacing-hPtop;

Prightstart = sz(1)-wPright-0.5*spacing; % Start of display to the right of the axes

hElement = 20*scalefact; % height of popup menu, checkboxes, small buttons
wButton1 = 60*scalefact;
hButton1 = 1.2*hElement;
wButton2 = wPright-spacing;
hButton2 = 80*scalefact;
hButton2b = 0.5*hButton2;
dh = 4*scalefact; % spacing (height)

ParTableh = hPtop-10;
ParTablex0 = spacing;
ParTabley0 = sz(2)-hPtop-spacing;
ParTablew = 765*scalefact+2*hElement+dh;
ParTableColw = (ParTablew-2*hElement+dh)/8;

Optionsx0 = ParTablex0+ParTablew+0.5*spacing;
Optionsy0 = ParTabley0+44*scalefact;
wOptionsLabel = 70*scalefact;
wOptionsSel = 145*scalefact;

Buttonsx0 = ParTablex0+ParTablew+wOptionsLabel+wOptionsSel+1.5*spacing;
Buttonsy0 = sz(2)-hPtop-spacing+dh;

Logx0 = Prightstart;
Logy0 = spacing;
Logw = wPright;
Logh = 110*scalefact;

FitSetx0 = Prightstart;
FitSety0 = Logy0+Logh+2*hElement;
FitSetw = wPright;
FitSeth = 125*scalefact;

Rmsdx0 = Prightstart;
Rmsdy0 = FitSety0+FitSeth+2*hElement;
Rmsdw = wPright;
Rmsdh = 125*scalefact;

% Axes
%-------------------------------------------------------------------------------
% Data display
hAx = axes('Parent',hFig,'Tag','dataaxes','Units','pixels',...
    'Position',[spacing spacing Axesw Axesh],'FontSize',8,'Layer','top');

NaNdata = NaN(1,numel(data));
mask = esfitdata.Opts.Mask;
dispData = esfitdata.data;
maxy = max(dispData(mask));
miny = min(dispData(mask));
YLimits = [miny maxy] + [-1 1]*Opt.PlotStretchFactor*(maxy-miny);
minx = min(esfitdata.Opts.x);
maxx = max(esfitdata.Opts.x);
x = esfitdata.Opts.x;

h(1) = line(x,NaNdata,'Color','k','Marker','.','LineStyle','none');
h(2) = line(x,NaNdata,'Color','r');
h(3) = line(x,NaNdata,'Color',[0 0.6 0]);
set(h(1),'Tag','expdata','XData',esfitdata.Opts.x,'YData',dispData);
set(h(2),'Tag','currsimdata');
set(h(3),'Tag','bestsimdata');
hAx.XLim = [minx maxx];
hAx.YLim = YLimits;
hAx.ButtonDownFcn = @axesButtonDownFcn;
grid(hAx,'on');
set(hAx,'YTickLabel',{})
if ~esfitdata.showxaxis
  set(hAx,'XTickLabel',{})
end
box on

showmaskedregions();

% Parameter table
%-------------------------------------------------------------------------------
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
uitable('Parent',hFig,'Tag','ParameterTable',...
    'FontSize',8,...
    'Position',[ParTablex0 ParTabley0 ParTablew ParTableh],...
    'ColumnFormat',columnformat,...
    'ColumnName',columnname,...
    'ColumnEditable',colEditable,...
    'CellEditCallback',@tableCellEditCallback,...
    'ColumnWidth',{hElement,hElement,ParTableColw,ParTableColw,ParTableColw,ParTableColw,ParTableColw,ParTableColw,ParTableColw,ParTableColw},...
    'RowName',[],...
    'Data',data,...
    'UserData',colEditable);
ParTableLabely0 = ParTabley0+ParTableh+dh;
uicontrol('Parent',hFig,'Style','text',...
    'Position',[ParTablex0 ParTableLabely0 2*wButton1 hElement],...
    'BackgroundColor',get(gcf,'Color'),...
    'FontWeight','bold','String','Parameters',...
    'HorizontalAl','left');

x0shift = ParTablew-4.5*wButton1-6*wButton1;
ParTableButtony0 = ParTableLabely0+dh/4;
uicontrol('Parent',hFig,'Style','text',...
    'Position',[ParTablex0+x0shift-0.2*wButton1 ParTableLabely0 1.2*wButton1 hElement],...
    'BackgroundColor',get(gcf,'Color'),...
    'FontWeight','bold','String','Start point:',...
    'HorizontalAl','left');
uicontrol('Parent',hFig,'Style','pushbutton','Tag','selectStartPointButtonCenter',...
    'Position',[ParTablex0+x0shift+wButton1 ParTableButtony0 wButton1 hButton1],...
    'String','center','Enable','on','Callback',@(src,evt) setStartPoint('center'),...
    'HorizontalAl','left',...
    'Tooltip','Set start values to center of range');
uicontrol('Parent',hFig,'Style','pushbutton','Tag','selectStartPointButtonRandom',...
    'Position',[ParTablex0+x0shift+2*wButton1 ParTableButtony0 wButton1 hButton1],...
    'String','random','Enable','on','Callback',@(src,evt) setStartPoint('random'),...
    'HorizontalAl','left',...
    'Tooltip','Set random start values');
uicontrol('Parent',hFig,'Style','pushbutton','Tag','selectStartPointButtonSelected',...
    'Position',[ParTablex0+x0shift+3*wButton1 ParTableButtony0 wButton1 hButton1],...
    'String','selected','Enable','on','Callback',@(src,evt) setStartPoint('selected'),...
    'HorizontalAl','left',...
    'Tooltip','Set start values from selected fit result');
uicontrol('Parent',hFig,'Style','pushbutton','Tag','selectStartPointButtonBest',...
    'Position',[ParTablex0+x0shift+4*wButton1 ParTableButtony0 wButton1 hButton1],...
    'String','best','Enable','on','Callback',@(src,evt) setStartPoint('best'),...
    'HorizontalAl','left',...
    'Tooltip','Set start values to current best fit');

  
x0shift = ParTablew-4*wButton1;
uicontrol('Parent',hFig,'Style','text',...
    'Position',[ParTablex0+x0shift-0.2*wButton1 ParTableLabely0 1.2*wButton1 hElement],...
    'BackgroundColor',get(gcf,'Color'),...
    'FontWeight','bold','String','Selection:',...
    'HorizontalAl','left');
uicontrol('Parent',hFig,'Style','pushbutton','Tag','selectInvButton',...
    'Position',[ParTablex0+x0shift+wButton1 ParTableButtony0 wButton1 hButton1],...
    'String','invert','Enable','on','Callback',@selectInvButtonCallback,...
    'HorizontalAl','left',...
    'Tooltip','Invert selection of parameters');
uicontrol('Parent',hFig,'Style','pushbutton','Tag','selectAllButton',...
    'Position',[ParTablex0+x0shift+2*wButton1 ParTableButtony0 wButton1 hButton1],...
    'String','all','Enable','on','Callback',@selectAllButtonCallback,...
    'HorizontalAl','left',...
    'Tooltip','Select all parameters');
uicontrol('Parent',hFig,'Style','pushbutton','Tag','selectNoneButton',...
    'Position',[ParTablex0+x0shift+3*wButton1 ParTableButtony0 wButton1 hButton1],...
    'String','none','Enable','on','Callback',@selectNoneButtonCallback,...
    'HorizontalAl','left',...
    'Tooltip','Unselect all parameters');

% FitOption selection
%-------------------------------------------------------------------------------
uicontrol('Parent',hFig,'Style','text',...
    'String','Function',...
    'Tooltip','Name of simulation function',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[Optionsx0 Optionsy0+5*(hElement+dh)+0.5*hElement wOptionsLabel hElement]);
uicontrol('Parent',hFig,'Style','text',...
    'String',esfitdata.fcnName,...
    'ForeGroundColor','b',...
    'Tooltip',sprintf('using output no. %d of %d',esfitdata.nOutArguments,esfitdata.OutArgument),...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[Optionsx0+wOptionsLabel Optionsy0+5*(hElement+dh)+0.5*hElement wOptionsSel hElement]);

uicontrol('Parent',hFig,'Style','text',...
    'String','Algorithm',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[Optionsx0 Optionsy0+4*(hElement+dh)-dh+0.5*hElement wOptionsLabel hElement]);
uicontrol('Parent',hFig,'Style','popupmenu',...
    'Tag','AlgorithmMenu',...
    'String',esfitdata.AlgorithmNames,...
    'Value',Opt.AlgorithmID,...
    'BackgroundColor','w',...
    'Tooltip','Fitting algorithm',...
    'Position',[Optionsx0+wOptionsLabel Optionsy0+4*(hElement+dh)+0.5*hElement wOptionsSel hElement]);

uicontrol('Parent',hFig,'Style','text',...
    'String','Target',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[Optionsx0 Optionsy0+3*(hElement+dh)-dh+0.5*hElement wOptionsLabel hElement]);
uicontrol('Parent',hFig,'Style','popupmenu',...
    'Tag','TargetMenu',...
    'String',esfitdata.TargetNames,...
    'Value',Opt.TargetID,...
    'BackgroundColor','w',...
    'Tooltip','Target function',...
    'Position',[Optionsx0+wOptionsLabel Optionsy0+3*(hElement+dh)+0.5*hElement wOptionsSel hElement]);

uicontrol('Parent',hFig,'Style','text',...
    'String','BaseLine',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[Optionsx0 Optionsy0+2*(dh+hElement)-dh+0.5*hElement wOptionsLabel hElement]);
uicontrol('Parent',hFig,'Style','popupmenu',...
    'Tag','BaseLineMenu',...
    'String',esfitdata.BaseLineStrings,...
    'Value',find(cellfun(@(x)x==esfitdata.BaseLine,esfitdata.BaseLineSettings),1),...
    'BackgroundColor','w',...
    'Tooltip','Baseline fitting',...
    'Position',[Optionsx0+wOptionsLabel Optionsy0+2*(dh+hElement)+0.5*hElement wOptionsSel hElement]);

uicontrol('Parent',hFig,'Style','text',...
    'String','AutoScale',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[Optionsx0 Optionsy0+(dh+hElement)-dh+0.5*hElement wOptionsLabel hElement]);
uicontrol('Parent',hFig,'Style','popupmenu',...
    'Tag','AutoScaleMenu',...
    'String',esfitdata.AutoScaleStrings,...
    'Value',find(cellfun(@(x)x==esfitdata.Opts.AutoScaleID,esfitdata.AutoScaleSettings),1),...
    'BackgroundColor','w',...
    'Tooltip','Autoscaling',...
    'Position',[Optionsx0+wOptionsLabel Optionsy0+(dh+hElement)+0.5*hElement wOptionsSel hElement]);

wMaskEl = 0.5*(wOptionsLabel+wOptionsSel);
uicontrol('Parent',hFig,'Style','checkbox',...
    'Tag','MaskCheckbox',...
    'String','Use mask',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'Value',1,...
    'BackgroundColor',get(gcf,'Color'),...
    'Tooltip','Use mask with excluded regions',...
    'Position',[Optionsx0+2*dh Optionsy0 wMaskEl hElement]);
uicontrol('Parent',hFig,'Style','pushbutton',...
    'Tag','SaveButton',...
    'String','Clear mask',...
    'Callback',@clearMaskCallback,...
    'Enable','on',...
    'Tooltip','Clear mask',...
    'Position',[Optionsx0+wMaskEl Optionsy0-dh wMaskEl hButton1]);

% Start/Stop buttons
%--------------------------------------------------------------------------
uicontrol('Parent',hFig,'Style','pushbutton',...
    'Tag','StartButton',...
    'String','Start fitting',...
    'Callback',@startButtonCallback,...
    'Visible','on',...
    'Tooltip','Start fitting',...
    'Position',[Buttonsx0 Buttonsy0-dh+3*hButton2b wButton2 hButton2]);
uicontrol('Parent',hFig,'Style','pushbutton',...
    'Tag','StopButton',...
    'String','Stop fitting',...
    'Visible','off',...
    'Tooltip','Stop fitting',...
    'Callback','global esfitdata; esfitdata.UserCommand = 1;',...
    'Position',[Buttonsx0 Buttonsy0-dh+3*hButton2b wButton2 hButton2]);
uicontrol('Parent',hFig,'Style','pushbutton',...
    'Tag','SaveButton',...
    'String','Save parameter set',...
    'Callback',@saveFitsetCallback,...
    'Enable','off',...
    'Tooltip','Save latest fitting result',...
    'Position',[Buttonsx0 Buttonsy0-dh+2*hButton2b wButton2 hButton2b]);
uicontrol('Parent',hFig,'Style','pushbutton',...
    'Tag','EvaluateButton',...
    'String','Evaluate at start point',...
    'Callback',@evaluateCallback,...
    'Enable','on',...
    'Tooltip','Run simulation for current start parameters',...
    'Position',[Buttonsx0 Buttonsy0-dh+hButton2b wButton2 hButton2b]);
uicontrol('Parent',hFig,'Style','pushbutton',...
    'Tag','ResetButton',...
    'String','Reset',...
    'Callback',@resetCallback,...
    'Enable','on',...
    'Tooltip','Clear fit history',...
    'Position',[Buttonsx0 Buttonsy0-dh wButton2 hButton2b]);

% Iteration and rmsd history displays
%-------------------------------------------------------------------------------
uicontrol('Parent',hFig,'Style','text',...
    'Position',[Rmsdx0 Rmsdy0+Rmsdh+4*hElement Rmsdw hElement],...
    'BackgroundColor',get(gcf,'Color'),...
    'FontWeight','bold','String','RMSD history',...
    'HorizontalAl','left');

h = uicontrol('Parent',hFig,'Style','text','Position',[Rmsdx0 Rmsdy0+Rmsdh+3*hElement 0.75*Rmsdw hElement]);
set(h,'FontSize',8,'String',' RMSD: -','ForegroundColor',[0 0.6 0],'Tooltip','Current best RMSD');
set(h,'Tag','RmsText','HorizontalAl','left');

h = uicontrol('Parent',hFig,'Style','checkbox','Position',[Rmsdx0+0.75*Rmsdw Rmsdy0+Rmsdh+3.2*hElement 0.30*Rmsdw hElement]);
set(h,'FontSize',8,'String','logscale','Tooltip','Set log scale on/off','Value',0);
set(h,'Tag','RmsLogPlot','Callback',@updatermsdplot);

hAx = axes('Parent',hFig,'Units','pixels','Position',[Rmsdx0 Rmsdy0+2.5*hElement Rmsdw-spacing Rmsdh],'Layer','top');
h = plot(hAx,1,NaN,'.');
set(h,'Tag','rmsdline','MarkerSize',5,'Color',[0.2 0.2 0.8]);
set(hAx,'FontSize',7,'YScale','lin','XTick',[],'YAxisLoc','right','Layer','top','YGrid','on');

h = uicontrol('Parent',hFig,'Style','text','Position',[Rmsdx0 Rmsdy0 0.9*Rmsdw 2*hElement]);
set(h,'FontSize',7,'Tag','logLine','Tooltip','Information from fitting algorithm');
set(h,'Horizontal','left');

% Fitset list
%-------------------------------------------------------------------------------
x0shift = 0;
wButton1 = FitSetw/4;
uicontrol('Parent',hFig,'Style','text','Tag','SetListTitle',...
    'Position',[FitSetx0 FitSety0+FitSeth+hElement FitSetw hElement],...
    'BackgroundColor',get(gcf,'Color'),...
    'FontWeight','bold','String','Parameter sets',...
    'Tooltip','List of stored fit parameter sets',...
    'HorizontalAl','left');
uicontrol('Parent',hFig,'Style','listbox','Tag','SetListBox',...
    'Position',[FitSetx0 FitSety0 FitSetw FitSeth],...
    'String','','Tooltip','',...
    'BackgroundColor',[1 1 0.9],...
    'KeyPressFcn',@deleteSetListKeyPressFcn,...
    'Callback',@setListCallback);
uicontrol('Parent',hFig,'Style','pushbutton','Tag','sortRMSDSetButton',...
    'Position',[FitSetx0+x0shift FitSety0+FitSeth+dh wButton1 hElement],...
    'String','rmsd',...
    'Tooltip','Sort parameter sets by rmsd','Enable','off',...
    'Callback',@sortRMSDSetButtonCallback);
uicontrol('Parent',hFig,'Style','pushbutton','Tag','sortIDSetButton',...
    'Position',[FitSetx0+x0shift+wButton1 FitSety0+FitSeth+dh wButton1 hElement],...
    'String','id',...
    'Tooltip','Sort parameter sets by ID','Enable','off',...
    'Callback',@sortIDSetButtonCallback);
uicontrol('Parent',hFig,'Style','pushbutton','Tag','exportSetButton',...
    'Position',[FitSetx0+x0shift+2*wButton1 FitSety0+FitSeth+dh wButton1 hElement],...
    'String','export',...
    'Tooltip','Export fit set to workspace','Enable','off',...
    'Callback',@exportSetButtonCallback);
uicontrol('Parent',hFig,'Style','pushbutton','Tag','deleteSetButton',...
    'Position',[FitSetx0+x0shift+3*wButton1 FitSety0+FitSeth+dh wButton1 hElement],...
    'String','delete',...
    'Tooltip','Delete fit set','Enable','off',...
    'Callback',@deleteSetButtonCallback);

% Error log panel
%-------------------------------------------------------------------------------
uicontrol('Parent',hFig,'Style','text',...
    'Position',[Logx0 Logy0+Logh Logw hElement],...
    'BackgroundColor',get(gcf,'Color'),...
    'FontWeight','bold','String','Log',...
    'Tooltip','Fitting information and error log',...
    'HorizontalAl','left');
hLogBox = uicontrol('Parent',hFig,'Style','listbox','Tag','LogBox',...
    'Position',[Logx0 Logy0 Logw Logh],...
    'String',{''},'Tooltip','',...
    'HorizontalAlignment','left',...
    'Min',0,'Max',2,...
    'Value',[],'Enable','inactive',...
    'BackgroundColor',[1 1 1]);

copymenu = uicontextmenu(hFig);

% Before R2017b (9.3), uimenu used Label instead of Text
if verLessThan('Matlab','9.3')
  menuTextProperty = 'Label';
else
  menuTextProperty = 'Text';
end
uimenu(copymenu,menuTextProperty,'Copy to clipboard','Callback',@copyLog);

% Before R2020a (9.8), uicontrol used UIContextMenu instead of ContextMenu
if verLessThan('Matlab','9.8')
  hLogBox.UIContextMenu = copymenu;
else
  hLogBox.ContextMenu = copymenu;
end

drawnow

set(hFig,'Visible','on')
set(hFig,'NextPlot','new');

end
%===============================================================================

%===============================================================================
function setStartPoint(sel)

global esfitdata

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
  case 'selected' % selected parameter set
    h = findobj('Tag','SetListBox');
    s = h.String;
    if ~isempty(s)
      s = s{h.Value};
      ID = sscanf(s,'%d');
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
hParamTable = findobj('Tag','ParameterTable');
data = get(hParamTable,'data');
for p = 1:numel(p_start)
  data{p,4} = sprintf('%0.6g',p_start(p));
  if updatebounds
    data{p,5} = sprintf('%0.6g',esfitdata.pvec_lb(p));
    data{p,6} = sprintf('%0.6g',esfitdata.pvec_ub(p));
  end
end
hParamTable.Data = data;

end
%===============================================================================

%===============================================================================
function axesButtonDownFcn(~,~)
global esfitdata
hAx = findobj('Tag','dataaxes');

% Get mouse-click point on axes
cp = hAx.CurrentPoint;
x = esfitdata.Opts.x;
x1 = cp(1,1);

% Create temporary patch updating with user mouse motion
maskColor = [1 1 1]*0.95;
tmpmask = patch(hAx,x1*ones(1,4),hAx.YLim([1 1 2 2]),maskColor,'Tag','maskPatch','EdgeColor','none');

% Move new patch to the back
c = hAx.Children([2:end 1]);
hAx.Children = c;

% Continuously update patch based on mouse position until next user click
set(gcf,'WindowButtonMotionFcn',@(hObject,eventdata) drawmaskedregion(tmpmask));
waitforbuttonpress;
set(gcf,'WindowButtonMotionFcn',[])

% Update masked regions
cp = hAx.CurrentPoint;
x2 = cp(1,1);
maskrange = sort([x1 x2]);
esfitdata.Opts.Mask(x>maskrange(1) & x<maskrange(2)) = 0;
delete(tmpmask);
showmaskedregions();
end
%===============================================================================

%===============================================================================
function drawmaskedregion(tmpmask)
cp = get (gca,'CurrentPoint');
xdata = tmpmask.XData;
xdata(2:3) = cp(1,1);
set(tmpmask,'XData',xdata);
end
%===============================================================================

%===============================================================================
function showmaskedregions()
global esfitdata
hAx = findobj('Tag','dataaxes');

% Delete existing mask patches
hMaskPatches = findobj(hAx,'Tag','maskPatch');
delete(hMaskPatches);

% Show masked-out regions
maskColor = [1 1 1]*0.95;
edges = find(diff([1; esfitdata.Opts.Mask(:); 1]));
upperlimit = numel(esfitdata.Opts.x);
edges(edges>upperlimit) = upperlimit;
excludedRegions = esfitdata.Opts.x(edges);
excludedRegions = reshape(excludedRegions,2,[]);

if ~isempty(excludedRegions)
  % Add a patch for each masked region
  nMaskPatches = size(excludedRegions,2);
  for r = 1:nMaskPatches
    x_patch = excludedRegions([1 2 2 1],r);
    y_patch = hAx.YLim([1 1 2 2]);
    patch(hAx,x_patch,y_patch,maskColor,'Tag','maskPatch','EdgeColor','none');
  end

  % Reorder so that mask patches are in the back
  c = hAx.Children([nMaskPatches+1:end, 1:nMaskPatches]);
  hAx.Children = c;
  drawnow
end

end
%===============================================================================
