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
%        .mask    array of 1 and 0 the same size as data vector
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

function result = esfit(data,fcn,p0,varargin)

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
    error('Opt.OutArg must contain two values [nOut iOut]');
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
    case 'iint',       Opt.TargetID = 3;
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
if ~isfield(Opt,'mask')
  Opt.mask = true(size(data_vec));
else
  if iscell(Opt.mask)
    Opt.mask = cat(2,Opt.mask{:});
  end
  Opt.mask = logical(Opt.mask(:));
  if numel(Opt.mask)~=numel(data_vec)
    error('Opt.mask has %d elements, but the data has %d elements.',numel(Opt.mask),numel(data));
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
  idx = esfitdata.idx{i};
  weights(idx) = Opt.weights(i);
end
esfitdata.weights = weights;

% x axis for plotting
if ~isfield(Opt,'x')
  Opt.x = 1:numel(esfitdata.data);
end

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
hFig = findall(0,'Tag','esfitFigure');
if ~isempty(hFig)
  close(hFig);
  clear hFig
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
      hPopup = findall(0,'Tag','algorithmpopup');
      set(findobj(hPopup,'Tag','nParamsField'),'Value',nActiveParams);
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
  end
else
  result.fit = fit;
  result.fitraw = fitraw;
end

result.baseline = baseline;
result.baselinetype = baselinetype;
result.scale = scale;
result.mask = esfitdata.Opts.mask;

result.bestfithistory.rmsd = esfitdata.besthistory.rmsd;
result.bestfithistory.pfit = esfitdata.besthistory.par;
if esfitdata.structureInputs
  result.bestfithistory.pfit2structs = esfitdata.p2args;
end

result.pnames = {esfitdata.pinfo.Name}.';
result.pnames = result.pnames(activeParams);
result.p_start = p_start;
result.pfit = pfit_active;
result.pfit_full = pfit;

result.pstd = pstd;
result.ci95 = ci95;
result.cov = covmatrix;
result.corr = corrmatrix;

result.argsfit = argsfit;

result.residuals = residuals0;
result.ssr = ssr0;
result.rmsd = rmsd0;
result.redchisquare = red_chisquare;

esfitdata.best = result;

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
  mask = Opt.mask;
else
  mask = true(size(Opt.mask));
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
    try
      % Get x-axis for EasySpin functions (and adapted functions based on EasySpin)
      [x,out{:}] = esfitdata.fcn(args{:});
      esfitdata.Opts.x = x;
    catch
      % Keep index axis for functions not returning an axis (e.g. curry)
      [out{:}] = esfitdata.fcn(args{:});
    end
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
rmsd = sqrt(mean(abs(residuals).^2.*esfitdata.weights(:)));
rmsd0 = sqrt(mean(abs(residuals0).^2));

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

if iscell(data)
  data = data{1};
end

global esfitdata
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
hFig = uifigure('Tag','esfitFigure');%,'WindowStyle','normal');
if ~strcmp(hFig.WindowStyle,'normal')
  hFig.WindowStyle = 'normal';
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
set(hFig,'MenuBar','none');
% set(hFig,'Resize','off');
set(hFig,'Name','esfit - Least-Squares Fitting','NumberTitle','off');
set(hFig,'CloseRequestFcn',@closeGUI);
  
spacing = 20*scalefact;
hPtop = 180*scalefact;
hElement = 22*scalefact; % height of popup menu, checkboxes, small buttons

fontsizetbl = 11;

% Set up grid layout
%-------------------------------------------------------------------------------
lgrid.main = uigridlayout(hFig,[1 2],'Padding',[1 1 1 1]*spacing,'ColumnSpacing',spacing);
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
hAx = uiaxes('Parent',lgrid.plot,'Tag','dataaxes','Layer','top');

NaNdata = NaN(1,numel(data));
mask = esfitdata.Opts.mask;
dispData = esfitdata.data;
maxy = max(dispData(mask));
miny = min(dispData(mask));
YLimits = [miny maxy] + [-1 1]*Opt.PlotStretchFactor*(maxy-miny);
minx = min(esfitdata.Opts.x);
maxx = max(esfitdata.Opts.x);
x = esfitdata.Opts.x;

h(1) = line(hAx,x,NaNdata,'Color',[1 1 1]*0.65);
h(2) = line(hAx,x,NaNdata,'Color','k','Marker','.','LineStyle','none');
h(3) = line(hAx,x,NaNdata,'Color','r');
h(4) = line(hAx,x,NaNdata,'Color',[0 0.6 0]);
set(h(1),'Tag','baselinedata');
set(h(1),'Visible','off');
set(h(2),'Tag','expdata','XData',esfitdata.Opts.x,'YData',dispData);
set(h(3),'Tag','currsimdata');
set(h(4),'Tag','bestsimdata');
hAx.XLim = [minx maxx];
hAx.YLim = YLimits;
hAx.ButtonDownFcn = @axesButtonDownFcn;
grid(hAx,'on');
box(hAx,'on')
set(hAx,'XTickLabel',{})

h(5) = yline(hAx,mean(dispData)-max(abs(dispData)));
set(h(5),'Tag','residualzero');
h(6) = line(hAx,x,NaNdata,'Color',[1 1 1]*0.65,'Marker','.','LineStyle','none');
set(h(6),'Tag','residualdata');
set(h(5),'Visible','off')
set(h(6),'Visible','off')

showmaskedregions();

% Fit parameter display
%-------------------------------------------------------------------------------
lgrid.parbanner = uigridlayout(lgrid.params,[1 4],'Padding',[0 0 0 0],'ColumnSpacing',20);
lgrid.parbanner.ColumnWidth = {'0.5x','1x','0.25x','1x'};
lgrid.parbanner.Layout.Row = 1;
uilabel('Parent',lgrid.parbanner,...
        'Text','Parameters',...
        'BackgroundColor',get(hFig,'Color'),...
        'FontWeight','bold',...
        'HorizontalAl','left','VerticalAl','center');

startpointbox = uigridlayout(lgrid.parbanner,[1 4],'Padding',[0 0 0 0],'ColumnSpacing',0);
uilabel('Parent',startpointbox,...
    'BackgroundColor',get(hFig,'Color'),...
    'FontWeight','bold','Text','Start point:',...
    'HorizontalAl','left','VerticalAl','center');
uibutton('Parent',startpointbox,'Tag','selectStartPointButtonCenter',...
    'Text','center','Enable','on','ButtonPushedFcn',@(src,evt) setStartPoint('center'),...
    'Tooltip','Set start values to center of range');
uibutton('Parent',startpointbox,'Tag','selectStartPointButtonRandom',...
    'Text','random','Enable','on','ButtonPushedFcn',@(src,evt) setStartPoint('random'),...
    'Tooltip','Set random start values');
uibutton('Parent',startpointbox,'Tag','selectStartPointButtonBest',...
    'Text','best','Enable','on','ButtonPushedFcn',@(src,evt) setStartPoint('best'),...
    'Tooltip','Set start values to current best fit');
  
selectbox = uigridlayout(lgrid.parbanner,[1 4],'Padding',[0 0 0 0],'ColumnSpacing',0);
selectbox.Layout.Column = 4;
uilabel('Parent',selectbox,...
    'BackgroundColor',get(hFig,'Color'),...
    'FontWeight','bold','Text','Selection:',...
    'HorizontalAl','left','VerticalAl','center');
uibutton('Parent',selectbox,'Tag','selectInvButton',...
    'Text','invert','Enable','on','ButtonPushedFcn',@(src,evt) selectButtonCallback('invert'),...
    'Tooltip','Invert selection of parameters');
uibutton('Parent',selectbox,'Tag','selectAllButton',...
    'Text','all','Enable','on','ButtonPushedFcn',@(src,evt) selectButtonCallback('all'),...
    'Tooltip','Select all parameters');
uibutton('Parent',selectbox,'Tag','selectNoneButton',...
    'Text','none','Enable','on','ButtonPushedFcn',@(src,evt) selectButtonCallback('none'),...
    'Tooltip','Unselect all parameters');

% Parameter table
columnname = {'','','Name','start','lower','upper','current','best','stdev','ci95 lower','ci95 upper'};
columnformat = {'char','logical','char','char','char','char','char','char','char','char','char'};
colEditable = [false true false true true true false false false false false];
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
  data{p,11} = '-';
end
hParamTable = uitable('Parent',lgrid.params,'Tag','ParameterTable',...
                      'ColumnFormat',columnformat,'ColumnName',columnname,'RowName',[],...
                      'ColumnEditable',colEditable,'FontSize',fontsizetbl,...
                      'CellEditCallback',@tableCellEditCallback,...
                      'ColumnWidth',{hElement,hElement,'1x','1x','1x','1x','1x','1x','1x','1x','1x'},...
                      'Data',data,...
                      'UserData',colEditable);
hParamTable.Layout.Row = 2;
if ~verLessThan('Matlab','9.7')
  s = uistyle('Interpreter','html');
  addStyle(hParamTable,s);
end

% FitOption selection
%-------------------------------------------------------------------------------
lgrid.fitcontrol = uigridlayout(lgrid.right,[1 2],'Padding',[0 0 0 0],'ColumnSpacing',spacing);
lgrid.fitcontrol.ColumnWidth = {'1.25x','1x'};

fitoptbox0 = uigridlayout(lgrid.fitcontrol,[3 1],'Padding',[0 0 0 0],'ColumnSpacing',0,'RowSpacing',5);
fitoptbox0.RowHeight = {hElement,'4x','2x'};

fitoptbox1 = uigridlayout(fitoptbox0,[1 2],'Padding',[0 0 0 0],'ColumnSpacing',0);
fitoptbox1.ColumnWidth = {'0.6x','0.8x',hElement};
uilabel('Parent',fitoptbox1,...
    'Text','Function',...
    'Tooltip','Name of simulation function',...
    'FontWeight','bold',...
    'HorizontalAl','left','VerticalAl','center',...
    'BackgroundColor',get(hFig,'Color'));
uilabel('Parent',fitoptbox1,...
    'Text',esfitdata.fcnName,...
    'FontColor','b',...
    'Tooltip',sprintf('using output no. %d of %d',esfitdata.nOutArguments,esfitdata.OutArgument),...
    'HorizontalAl','left','VerticalAl','center',...
    'BackgroundColor',get(hFig,'Color'));
uibutton('Parent',fitoptbox1,'Tag','AlgorithmSettingsButton',...
         'Text','','Tooltip','Set fitting algorithm options',...
         'Icon','private/settingsicon.png','IconAlignment','center',...
         'Enable','on');

fitoptbox2 = uigridlayout(fitoptbox0,[4 2],'Padding',[0 0 0 0],'ColumnSpacing',0,'RowSpacing',5);
fitoptbox2.ColumnWidth = {'0.6x','1x'};
uilabel('Parent',fitoptbox2,...
    'Text','Algorithm',...
    'FontWeight','bold',...
    'HorizontalAlign','left','VerticalAlign','center',...
    'BackgroundColor',get(hFig,'Color'));
uidropdown('Parent',fitoptbox2,...
    'Tag','AlgorithmMenu',...
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
    'BackgroundColor',get(hFig,'Color'));
uidropdown('Parent',fitoptbox2,...
    'Tag','TargetMenu',...
    'Items',esfitdata.TargetNames,...
    'ItemsData',1:numel(esfitdata.TargetNames),...
    'Value',Opt.TargetID,...
    'BackgroundColor','w',...
    'Tooltip','Target function');

uilabel('Parent',fitoptbox2,...
    'Text','AutoScale',...
    'FontWeight','bold',...
    'HorizontalAlign','left','VerticalAlign','center',...
    'BackgroundColor',get(hFig,'Color'));
uidropdown('Parent',fitoptbox2,...
    'Tag','AutoScaleMenu',...
    'Items',esfitdata.AutoScaleStrings,...
    'ItemsData',esfitdata.AutoScaleSettings,...
    'Value',esfitdata.Opts.AutoScaleID,...
    'BackgroundColor','w',...
    'Tooltip','Autoscaling');

uilabel('Parent',fitoptbox2,...
    'Text','BaseLine',...
    'FontWeight','bold',...
    'HorizontalAlign','left','VerticalAlign','center',...
    'BackgroundColor',get(hFig,'Color'));
uidropdown('Parent',fitoptbox2,...
    'Tag','BaseLineMenu',...
    'Items',esfitdata.BaseLineStrings,...
    'ItemsData',esfitdata.BaseLineSettings,...
    'Value',esfitdata.BaseLine,...
    'BackgroundColor','w',...
    'Tooltip','Baseline fitting');

fitoptbox3 = uigridlayout(fitoptbox0,[2 2],'Padding',[0 0 0 0],'ColumnSpacing',5,'RowSpacing',5);
uicheckbox('Parent',fitoptbox3,'Tag','BaselineCheckbox',...
           'Text','plot baseline','Tooltip','Plot baseline on top of data',...
           'Value',0,'ValueChangedFcn',@showbaseline);
uicheckbox('Parent',fitoptbox3,'Tag','ResidualCheckbox',...
           'Text','plot residuals','Tooltip','Plot residuals',...
           'Value',0,'ValueChangedFcn',@showresiduals);

uicheckbox('Parent',fitoptbox3,'Tag','MaskCheckbox',...
           'Text','use mask','Tooltip','Use mask with excluded regions',...
           'Value',1);
uibutton('Parent',fitoptbox3,'Tag','SaveButton',...
         'Text','clear mask','Tooltip','Clear mask',...
         'Enable','on',...
         'ButtonPushedFcn',@clearMaskCallback);

% Start/Stop buttons
%--------------------------------------------------------------------------
fitctrlbox = uigridlayout(lgrid.fitcontrol,[4 1],'Padding',[0 0 0 0],'RowSpacing',0);
fitctrlbox.RowHeight = {'2x','1x','1x','1x'};
uibutton('Parent',fitctrlbox,...
    'Tag','StartButton',...
    'Text','Start fitting',...
    'ButtonPushedFcn',@startButtonCallback,...
    'Visible','on',...
    'Tooltip','Start fitting');
uibutton('Parent',fitctrlbox,...
    'Tag','SaveButton',...
    'Text','Save parameter set',...
    'ButtonPushedFcn',@saveFitsetCallback,...
    'Enable','off',...
    'Tooltip','Save latest fitting result');
uibutton('Parent',fitctrlbox,...
    'Tag','EvaluateButton',...
    'Text','Evaluate at start point',...
    'ButtonPushedFcn',@evaluateCallback,...
    'Enable','on',...
    'Tooltip','Run simulation for current start parameters');
uibutton('Parent',fitctrlbox,...
    'Tag','ResetButton',...
    'Text','Reset',...
    'ButtonPushedFcn',@resetCallback,...
    'Enable','on',...
    'Tooltip','Clear fit history');

% Fitset list
%-------------------------------------------------------------------------------
lgrid.fitsets = uigridlayout(lgrid.right,[2 1],'Padding',[0 0 0 0],'RowSpacing',5);
lgrid.fitsets.RowHeight = {hElement,'1x',hElement};

uilabel('Parent',lgrid.fitsets,'Tag','SetListTitle',...
    'BackgroundColor',get(hFig,'Color'),...
    'FontWeight','bold','Text','Fit parameter sets',...
    'Tooltip','List of stored fit parameter sets',...
    'HorizontalAl','left','VerticalAl','top');

% Fit parameter set table
columnname = {'ID','rmsd','scale','baseline',''};
columnformat = {'char','char','char','char','char'};
colEditable = [false false false false false];
hFitSetTable = uitable('Parent',lgrid.fitsets,'Tag','FitSetTable',...
                      'ColumnFormat',columnformat,'ColumnName',columnname,'RowName',[],...
                      'ColumnEditable',colEditable,'FontSize',fontsizetbl,...
                      'ColumnWidth',{1.5*hElement,'1x','1x','1x','1x'},...
                      'Data',[],...
                      'SelectionChangedFcn',@setListCallback,...
                      'Enable','on',...
                      'KeyPressFcn',@deleteSetListKeyPressFcn);
if ~verLessThan('matlab','9.11')  % 9.11 = R2021b
  hFitSetTable.SelectionType = 'row';
end
if ~verLessThan('matlab','9.7')  % 9.7 = R2019b
  hFitSetTable.ColumnSortable = true;
end
fitsetbuttonbox = uigridlayout(lgrid.fitsets,[1 5],'Padding',[0 0 0 0],'ColumnSpacing',0);
fitsetbuttonbox.ColumnWidth = {2*hElement,'1x','1x','1x',2*hElement};
h(1) = uibutton('Parent',fitsetbuttonbox,'Tag','selectStartPointButtonSelected',...
    'Text','set as start',...
    'Tooltip','Set as start point for fitting','Enable','off',...
    'ButtonPushedFcn',@(src,evt) setStartPoint('selected'));
h(1).Layout.Column = 2;
h(2) = uibutton('Parent',fitsetbuttonbox,'Tag','exportSetButton',...
    'Text','export',...
    'Tooltip','Export fit set to workspace','Enable','off',...
    'ButtonPushedFcn',@exportSetButtonCallback);
h(2).Layout.Column = 3;
h(3) = uibutton('Parent',fitsetbuttonbox,'Tag','deleteSetButton',...
    'Text','delete',...
    'Tooltip','Delete fit set','Enable','off',...
    'ButtonPushedFcn',@deleteSetButtonCallback);
h(3).Layout.Column = 4;

% Iteration and rmsd history displays
%-------------------------------------------------------------------------------
lgrid.rmsdplot = uigridlayout(lgrid.right,[4 1],'Padding',[0 0 0 0],'RowSpacing',0);
lgrid.rmsdplot.RowHeight = {hElement,hElement,'1x',hElement};

uilabel('Parent',lgrid.rmsdplot,...
    'BackgroundColor',get(hFig,'Color'),...
    'FontWeight','bold','Text','RMSD history',...
    'HorizontalAl','left','VerticalAl','center');

rmsdinfobox = uigridlayout(lgrid.rmsdplot,[1 2],'Padding',[0 0 0 0]);
rmsdinfobox.ColumnWidth = {'1x',hElement};
uilabel('Parent',rmsdinfobox,'Tag','RmsText',...
        'Text',' RMSD: -',...
        'FontColor',[0 0.6 0],'Tooltip','Current best RMSD',...
        'HorizontalAl','left','VerticalAl','center');

uibutton('state','Parent',rmsdinfobox,'Tag','RmsLogPlot',...
         'Text','log','FontSize',fontsizetbl-1,'Tooltip','Set log scale on/off',...
         'Value',0,'ValueChangedFcn',@updatermsdplot);

hAx = uiaxes('Parent',lgrid.rmsdplot,'Units','pixels','Layer','top');
h = plot(hAx,1,NaN,'.');
set(h,'Tag','rmsdline','MarkerSize',5,'Color',[0.2 0.2 0.8]);
set(hAx,'box','on','YScale','lin','XTick',[],'YAxisLoc','right','Layer','top','YGrid','on');

uilabel('Parent',lgrid.rmsdplot,'Tag','logLine',...
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
    'BackgroundColor',get(hFig,'Color'),...
    'FontWeight','bold','Text','Log',...
    'Tooltip','Fitting information and error log',...
    'HorizontalAl','left','VerticalAl','top');
uibutton('Parent',logtitlebox,'Tag','copyLogButton',...
         'Text','','Tooltip','Copy to clipboard',...
         'Icon','private/copyicon.png','IconAlignment','center',...
         'Enable','on',...
         'ButtonPushedFcn',@copyLog);

hLogBox = uitextarea('Parent',lgrid.log,'Tag','LogBox',...
    'Value',{''},'Tooltip','',...
    'FontSize',fontsizetbl,...
    'WordWrap','on','Editable','off',...
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

% Algorithm settings pop-up figure
%-------------------------------------------------------------------------------
hPopup = uifigure('Tag','algorithmpopup','WindowStyle','normal');
set(hPopup,'Visible','off')

sz = [400 320]; % popup size
sz = sz*scalefact;
xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the screen vertically
set(hPopup,'position',[xpos, ypos, sz(1), sz(2)],'units','pixels');
set(hPopup,'MenuBar','none');
% set(hFig,'Resize','off');
set(hPopup,'Name','esfit - Algorithm settings','NumberTitle','off');
set(hPopup,'CloseRequestFcn',...
    'hPopup = findall(0,''Tag'',''algorithmpopup''); set(hPopup,''Visible'',''off'');');
  
% Set up drop down menu and tabs
%-------------------------------------------------------------------------------
lgrid.main = uigridlayout(hPopup,[3 1],'Padding',[1 1 1 1]*spacing,'RowSpacing',5);
lgrid.main.RowHeight = {hElement,hElement,'1x'};

dropdownbox = uigridlayout(lgrid.main,[1 4],'Padding',[0 0 0 0],'ColumnSpacing',0);
dropdownbox.ColumnWidth = {'0.6x','1x','0.2x','0.6x'};
uilabel('Parent',dropdownbox,...
    'Text','Algorithm',...
    'FontWeight','bold',...
    'HorizontalAlign','left','VerticalAlign','center',...
    'BackgroundColor',get(hPopup,'Color'));
uidropdown('Parent',dropdownbox,...
    'Tag','AlgorithmMenuPopup',...
    'Items',esfitdata.AlgorithmNames,...
    'ItemsData',1:numel(esfitdata.AlgorithmNames),...
    'Value',Opt.AlgorithmID,...
    'BackgroundColor','w',...
    'ValueChangedFcn',@(src,evt) selectAlgorithm(src),...
    'Tooltip','Fitting algorithm');
b = uibutton('Parent',dropdownbox,...
         'Tag','setAlgorithmDefaults',...
         'Text','Set to defaults','Enable','on','ButtonPushedFcn',@(src,evt) updateAlgorithmDefaults,...
         'Tooltip','Reset all parameters to the default');
b.Layout.Column = 4;

nparamsbox = uigridlayout(lgrid.main,[1 4],'Padding',[0 0 0 0],'ColumnSpacing',0);
nparamsbox.ColumnWidth = {'1.2x','0.4x','0.4x','1x'};
l = uilabel('Parent',nparamsbox,...
    'Text','N',...
    'FontWeight','bold',...
    'HorizontalAlign','center','VerticalAlign','center',...
    'BackgroundColor',get(hPopup,'Color'));
l.Layout.Column = 2;
uieditfield('numeric','Parent',nparamsbox,...
        'Tag','nParamsField',...
        'ToolTip','number of fitting parameters',...
        'Value',esfitdata.nParameters,...
        'Editable','off');

tabs = uitabgroup(lgrid.main,'TabLocation','left','Tag','AlgorithmTabs');
AlgorithmAbbrev{1} = 'simplex';
AlgorithmAbbrev{2} = 'levmar';
AlgorithmAbbrev{3} = 'montecarlo';
AlgorithmAbbrev{4} = 'genetic';
AlgorithmAbbrev{5} = 'grid';
AlgorithmAbbrev{6} = 'swarm';
for i = 1:(numel(esfitdata.AlgorithmNames)-1)
  t(i) = uitab(tabs,'Title',AlgorithmAbbrev{i},'ToolTip',esfitdata.AlgorithmNames{i});
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
        'BackgroundColor',get(hPopup,'Color'));
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
  tabs.UserData{i} = hsetting;
end
updateAlgorithmDefaults()

% Set callback for settings button to open popup menu
settingsbutton = findobj(hFig,'Tag','AlgorithmSettingsButton');
set(settingsbutton,'ButtonPushedFcn',@openAlgorithmSettings);

set(hFig,'Visible','on')
set(hFig,'NextPlot','new');
set(hPopup,'NextPlot','new');

end
%===============================================================================

%===============================================================================
function evaluateCallback(~,~)
% Evaluate for selected parameters
global esfitdata
hFig = findall(0,'Tag','esfitFigure');

esfitdata.modelErrorHandler = @(ME) GUIErrorHandler(ME);

p_eval = esfitdata.p_start;
active = ~esfitdata.fixedParams;
p_eval = p_eval(active);
expdata = esfitdata.data(:);
esfitdata.Opts.AutoScaleID = get(findobj(hFig,'Tag','AutoScaleMenu'),'Value');
esfitdata.Opts.BaseLine = get(findobj(hFig,'Tag','BaseLineMenu'),'Value');
esfitdata.Opts.useMask = get(findobj(hFig,'Tag','MaskCheckbox'),'Value')==1;
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
currbaseline = real(esfitdata.curr.baseline(:));

% Update plotted data
x = esfitdata.Opts.x(:);
set(findobj(hFig,'Tag','expdata'),'XData',x,'YData',expdata);
set(findobj(hFig,'Tag','currsimdata'),'XData',x,'YData',currsim);
set(findobj(hFig,'Tag','baselinedata'),'XData',x,'YData',currbaseline);
residualzero = get(findobj(hFig,'Tag','residualzero'),'Value');
residualdata = (expdata-currsim)+residualzero;
set(findobj(hFig,'Tag','residualdata'),'XData',x,'YData',residualdata);

% Readjust vertical range
updateaxislimits()
drawnow

% Readjust mask patches
maskPatches = findobj(hFig,'Tag','maskPatch');
for mp = 1:numel(maskPatches)
  maskPatches(mp).YData = YLimits([1 1 2 2]).';
end

% Update column with best values if current parameter set is new best
str = sprintf('Current RMSD: %g\n',rmsd);
hRmsText = findobj(hFig,'Tag','RmsText');
set(hRmsText,'Text',str,'FontColor',[1 0 0]);

end
%===============================================================================

%===============================================================================
function startButtonCallback(~,~)

global esfitdata
hFig = findall(0,'Tag','esfitFigure');
hPopup = findall(0,'Tag','algorithmpopup');

switch get(findobj(hFig,'Tag','StartButton'),'Text')
  case 'Start fitting'
    esfitdata.UserCommand = 0;
  case 'Stop fitting'
    esfitdata.UserCommand = 1;
    return
end

% Update GUI
%-------------------------------------------------------------------------------

% Hide Start button, show Stop button
set(findobj(hFig,'Tag','StartButton'),'Text','Stop fitting');
set(findobj(hFig,'Tag','SaveButton'),'Enable','off');

% Disable other buttons
set(findobj(hFig,'Tag','EvaluateButton'),'Enable','off');
set(findobj(hFig,'Tag','ResetButton'),'Enable','off');
set(findobj(hFig,'Tag','AlgorithmSettingsButton'),'Enable','off');

% Disable listboxes
set(findobj(hFig,'Tag','AlgorithmMenu'),'Enable','off');
set(findobj(hFig,'Tag','TargetMenu'),'Enable','off');
set(findobj(hFig,'Tag','AutoScaleMenu'),'Enable','off');
set(findobj(hFig,'Tag','BaseLineMenu'),'Enable','off');

% Disable parameter table
set(findobj(hFig,'Tag','selectAllButton'),'Enable','off');
set(findobj(hFig,'Tag','selectNoneButton'),'Enable','off');
set(findobj(hFig,'Tag','selectInvButton'),'Enable','off');
set(findobj(hFig,'Tag','selectStartPointButtonCenter'),'Enable','off');
set(findobj(hFig,'Tag','selectStartPointButtonRandom'),'Enable','off');
set(findobj(hFig,'Tag','selectStartPointButtonBest'),'Enable','off');
colEditable = get(findobj(hFig,'Tag','ParameterTable'),'UserData');
set(findobj(hFig,'Tag','ParameterTable'),'ColumnEditable',false(size(colEditable)));
set(findobj(hFig,'Tag','ParameterTable'),'CellEditCallback',[]);

% Disable algorithm settings window
set(findobj(hPopup,'Tag','AlgorithmMenuPopup'),'Enable','off');
set(findobj(hPopup,'Tag','setAlgorithmDefaults'),'Enable','off');

% Remove displayed best fit and uncertainties
hTable = findobj(hFig,'Tag','ParameterTable');
Data = hTable.Data;
for p = 1:size(Data,1)
  Data{p,7} = '-';
  Data{p,8} = '-';
  Data{p,9} = '-';
  Data{p,10} = '-';
  Data{p,11} = '-';
end
set(hTable,'Data',Data);

% Get fixed parameters
for p = 1:esfitdata.nParameters
  esfitdata.fixedParams(p) = Data{p,2}==0;
end

% Disable fitset list controls
set(findobj(hFig,'Tag','selectStartPointButtonSelected'),'Enable','off');
set(findobj(hFig,'Tag','exportSetButton'),'Enable','off');
set(findobj(hFig,'Tag','deleteSetButton'),'Enable','off');

% Disable mask tools
hAx = findobj(hFig,'Tag','dataaxes');
hAx.ButtonDownFcn = [];
set(findobj(hFig,'Tag','clearMaskButton'),'Enable','off');
set(findobj(hFig,'Tag','MaskCheckbox'),'Enable','off');

% Pull settings from UI
%-------------------------------------------------------------------------------
% Determine selected method, target, autoscaling, start point
esfitdata.Opts.AlgorithmID = get(findobj(hFig,'Tag','AlgorithmMenu'),'Value');
esfitdata.Opts.TargetID = get(findobj(hFig,'Tag','TargetMenu'),'Value');
esfitdata.Opts.AutoScaleID = get(findobj(hFig,'Tag','AutoScaleMenu'),'Value');
esfitdata.Opts.BaseLine = get(findobj(hFig,'Tag','BaseLineMenu'),'Value');
esfitdata.Opts.useMask = get(findobj(hFig,'Tag','MaskCheckbox'),'Value')==1;

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
esfitdata.currFitSet.Mask = esfitdata.Opts.useMask && ~all(esfitdata.Opts.mask);


% Update GUI with fit results
%-------------------------------------------------------------------------------

% Remove current values and uncertainties from parameter table
hTable = findobj(hFig,'Tag','ParameterTable');
Data = hTable.Data;
pi = 1;
for p = 1:size(Data,1)
  Data{p,7} = '-'; 
  if ~esfitdata.fixedParams(p) && ~isempty(esfitdata.best.pstd)
    Data{p,8} = sprintf('%0.6g',esfitdata.best.pfit(pi));
    Data{p,9} = sprintf('%0.6g',esfitdata.best.pstd(pi));
    Data{p,10} = sprintf('%0.6g',esfitdata.best.ci95(pi,1));
    Data{p,11} = sprintf('%0.6g',esfitdata.best.ci95(pi,2));
    pi = pi+1;
  else
    Data{p,9} = '-';
    Data{p,10} = '-';
    Data{p,11} = '-';
  end
end
set(hTable,'Data',Data);

% Hide current sim plot in data axes
set(findobj(hFig,'Tag','currsimdata'),'YData',NaN(1,numel(esfitdata.data)));
residualzero = get(findobj(hFig,'Tag','residualzero'),'Value');
set(findobj(hFig,'Tag','residualdata'),'YData',(esfitdata.data-esfitdata.best.fit)+residualzero);
drawnow

% Reactivate UI components
set(findobj(hFig,'Tag','SaveButton'),'Enable','on');

if isfield(esfitdata,'FitSets') && numel(esfitdata.FitSets)>0
  set(findobj(hFig,'Tag','selectStartPointButtonSelected'),'Enable','on');
  set(findobj(hFig,'Tag','exportSetButton'),'Enable','on');
  set(findobj(hFig,'Tag','deleteSetButton'),'Enable','on');
end

% Hide stop button, show start button
set(findobj(hFig,'Tag','StartButton'),'Text','Start fitting');

% Re-enable other buttons
set(findobj(hFig,'Tag','EvaluateButton'),'Enable','on');
set(findobj(hFig,'Tag','ResetButton'),'Enable','on');
set(findobj(hFig,'Tag','AlgorithmSettingsButton'),'Enable','on');

% Re-enable listboxes
set(findobj(hFig,'Tag','AlgorithmMenu'),'Enable','on');
set(findobj(hFig,'Tag','TargetMenu'),'Enable','on');
set(findobj(hFig,'Tag','AutoScaleMenu'),'Enable','on');
set(findobj(hFig,'Tag','BaseLineMenu'),'Enable','on');
set(findobj(hFig,'Tag','AutoScaleMenu'),'Enable','on');

% Re-enable parameter table and its selection controls
set(findobj(hFig,'Tag','selectAllButton'),'Enable','on');
set(findobj(hFig,'Tag','selectNoneButton'),'Enable','on');
set(findobj(hFig,'Tag','selectInvButton'),'Enable','on');
set(findobj(hFig,'Tag','selectStartPointButtonCenter'),'Enable','on');
set(findobj(hFig,'Tag','selectStartPointButtonRandom'),'Enable','on');
set(findobj(hFig,'Tag','selectStartPointButtonBest'),'Enable','on');
set(findobj(hFig,'Tag','ParameterTable'),'ColumnEditable',colEditable);
set(findobj(hFig,'Tag','ParameterTable'),'CellEditCallback',@tableCellEditCallback);

% Re-enable mask tools
hAx = findobj(hFig,'Tag','dataaxes');
hAx.ButtonDownFcn = @axesButtonDownFcn;
set(findobj(hFig,'Tag','clearMaskButton'),'Enable','on');
set(findobj(hFig,'Tag','MaskCheckbox'),'Enable','on');

% Re-enable algorithm settings window
set(findobj(hPopup,'Tag','AlgorithmMenuPopup'),'Enable','on');
set(findobj(hPopup,'Tag','setAlgorithmDefaults'),'Enable','on');

end
%===============================================================================

%===============================================================================
function resetCallback(~,~)
% Reset best fit
global esfitdata
hFig = findall(0,'Tag','esfitFigure');

% Remove messages from log
scroll(findobj(hFig,'Tag','LogBox'),'top')
set(findobj(hFig,'Tag','LogBox'),'Value',{''})

% Remove best fit simulation from plot
hBestSim = findobj(hFig,'Tag','bestsimdata');
hBestSim.YData = NaN(size(hBestSim.YData));
esfitdata.best = [];

% Remove baseline from plot
hBaseLine = findobj(hFig,'Tag','baselinedata');
hBaseLine.YData = NaN(size(hBaseLine.YData));

% Remove current fit simulation from plot
hCurrSim = findobj(hFig,'Tag','currsimdata');
hCurrSim.YData = NaN(size(hCurrSim.YData));
esfitdata.curr = [];

% Remove residuals from plot
hResiduals = findobj(hFig,'Tag','residualdata');
hResiduals.YData = NaN(size(hResiduals.YData));

% Readjust vertical range
mask = esfitdata.Opts.mask;
expdata = esfitdata.data(:);
maxy = max(expdata(mask));
miny = min(expdata(mask));
YLimits = [miny maxy] + [-1 1]*esfitdata.Opts.PlotStretchFactor*(maxy-miny);
set(findobj(hFig,'Tag','dataaxes'),'YLim',YLimits);
drawnow

% Reset rmsdhistory plot
esfitdata.rmsdhistory = [];
updatermsdplot;
iterationprint('');

% Reset besthistory 
esfitdata.besthistory.rmsd = [];
esfitdata.besthistory.par = [];

% Remove displayed best fit and uncertainties
hTable = findobj(hFig,'Tag','ParameterTable');
Data = hTable.Data;
for p = 1:size(Data,1)
  Data{p,7} = '-';
  Data{p,8} = '-';
  Data{p,9} = '-';
  Data{p,10} = '-';
  Data{p,11} = '-';
end
set(hTable,'Data',Data);

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

hFig = findall(0,'Tag','esfitFigure');

% Get relevant quantities
x = esfitdata.Opts.x(:);
expdata = esfitdata.data(:);
bestsim = real(esfitdata.best.fit(:));
currsim = real(esfitdata.curr.sim(:));
currbaseline = real(esfitdata.curr.baseline(:));
currpar = esfitdata.curr.par;
bestpar = esfitdata.best.par;

% Update plotted data
set(findobj(hFig,'Tag','expdata'),'XData',x,'YData',expdata);
set(findobj(hFig,'Tag','bestsimdata'),'XData',x,'YData',bestsim);
set(findobj(hFig,'Tag','currsimdata'),'XData',x,'YData',currsim);
set(findobj(hFig,'Tag','baselinedata'),'XData',x,'YData',currbaseline);
residualzero = get(findobj(hFig,'Tag','residualzero'),'Value');
residualdata = (expdata-currsim)+residualzero;
set(findobj(hFig,'Tag','residualdata'),'XData',x,'YData',residualdata);

% Readjust vertical range
updateaxislimits()

% Readjust mask patches
maskPatches = findobj(hFig,'Tag','maskPatch');
for mp = 1:numel(maskPatches)
  maskPatches(mp).YData = YLimits([1 1 2 2]).';
end

% Update column with current parameter values
hParamTable = findobj(hFig,'Tag','ParameterTable');
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
  hRmsText = findobj(hFig,'Tag','RmsText');
  set(hRmsText,'Text',str,'FontColor',[0 0.6 0]);

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
hFig = findall(0,'Tag','esfitFigure');
% Update rmsd plot
hRmsText = findobj(hFig,'Tag','RmsText');
if isfield(esfitdata,'best') && ~isempty(esfitdata.best) && ~isempty(esfitdata.best.rmsd)
  str = sprintf('Current best RMSD: %g\n',esfitdata.best.rmsd);
else
  str = sprintf('Current best RMSD: -\n');
end
set(hRmsText,'Text',str);

hRmsLogPlot = findobj(hFig,'Tag','RmsLogPlot');

hrmsdline = findobj(hFig,'Tag','rmsdline');
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
function iterationprint(str)
hFig = findall(0,'Tag','esfitFigure');
hLogLine = findobj(hFig,'Tag','logLine');
if isempty(hLogLine)
  disp(strtrim(str));
else
  set(hLogLine,'Text',strtrim(str));
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
hFig = findall(0,'Tag','esfitFigure');
h = findobj(hFig,'Tag','FitSetTable');
if isempty(h.Selection), return; end
data = h.Data;
if ~isempty(data)
  k = find([esfitdata.FitSets.ID]==data{h.Selection,1});
  if k>0
    fitset = esfitdata.FitSets(k);

    % Set column with best-fit parameter values
    hTable = findobj(hFig,'Tag','ParameterTable');
    data = get(hTable,'data');

    pi = 1;
    for p = 1:size(data,1)
      data{p,8} = sprintf('%0.6g',fitset.pfit_full(p));
      if ~fitset.fixedParams(p) && ~isempty(fitset.pstd)
        data{p,9} = sprintf('%0.6g',fitset.pstd(pi));
        data{p,10} = sprintf('%0.6g',fitset.ci95(pi,1));
        data{p,11} = sprintf('%0.6g',fitset.ci95(pi,2));
        pi = pi+1;
      else
        data{p,9} = '-';
        data{p,10} = '-';
        data{p,11} = '-';
      end
    end
    set(hTable,'Data',data);

    h = findobj(hFig,'Tag','bestsimdata');
    set(h,'YData',fitset.fit);
    h = findobj(hFig,'Tag','baselinedata');
    set(h,'YData',fitset.baseline);
    drawnow
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
function refreshFitsetList(idx)
global esfitdata
hFig = findall(0,'Tag','esfitFigure');
h = findobj(hFig,'Tag','FitSetTable');
nSets = numel(esfitdata.FitSets);
data = cell(nSets,4);
for k = 1:nSets
  data{k,1} = esfitdata.FitSets(k).ID;
  data{k,2} = esfitdata.FitSets(k).rmsd;
  data{k,3} = esfitdata.FitSets(k).scale(1);
  data{k,4} = esfitdata.FitSets(k).baselinetype{1};
  if esfitdata.FitSets(k).Mask
    maskstr = ' (mask)';
  else
    maskstr = '';
  end
  data{k,5} = maskstr;
end
set(h,'Data',data);
if idx>0, set(h,'Selection',idx); end
if idx==-1, set(h,'Selection',nSets); end

if nSets>0, state = 'on'; else, state = 'off'; end
set(findobj(hFig,'Tag','selectStartPointButtonSelected'),'Enable',state);
set(findobj(hFig,'Tag','exportSetButton'),'Enable',state);
set(findobj(hFig,'Tag','deleteSetButton'),'Enable',state);

displayFitSet;
end
%===============================================================================

%===============================================================================
function exportSetButtonCallback(~,~)
global esfitdata
hFig = findall(0,'Tag','esfitFigure');
h = findobj(hFig,'Tag','FitSetTable');
ID = h.Data{h.Selection,1};
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
function deleteSetButtonCallback(~,~)
global esfitdata
hFig = findall(0,'Tag','esfitFigure');
h = findobj(hFig,'Tag','FitSetTable');
idx = h.Selection;
nSets = numel(esfitdata.FitSets);
if nSets>0
    ID = h.Data{h.Selection,1};
    for k = nSets:-1:1
        if esfitdata.FitSets(k).ID==ID
            esfitdata.FitSets(k) = [];
        end
    end
    if idx>length(esfitdata.FitSets), idx = length(esfitdata.FitSets); end
    if idx==0, idx = 1; end
    h.Selection = idx;
    refreshFitsetList(0);
end

if isempty(h.Data)
    set(findobj(hFig,'Tag','selectStartPointButtonSelected'),'Enable','off');
    set(findobj(hFig,'Tag','exportSetButton'),'Enable','off');
    set(findobj(hFig,'Tag','deleteSetButton'),'Enable','off');
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
function setStartPoint(sel)

global esfitdata
hFig = findall(0,'Tag','esfitFigure');

% Update fixed parameters
hTable = findobj(hFig,'Tag','ParameterTable');
Data = hTable.Data;
for p = 1:esfitdata.nParameters
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
    h = findobj(hFig,'Tag','FitSetTable');
    data = h.Data;
    if ~isempty(data)
      ID = data{h.Selection,1};
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
hParamTable = findobj(hFig,'Tag','ParameterTable');
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
function selectButtonCallback(type)
hFig = findall(0,'Tag','esfitFigure');
h = findobj(hFig,'Tag','ParameterTable');
d = h.Data;
switch type
  case 'all'
    d(:,2) = {true};
  case 'none'
    d(:,2) = {false};
  case 'invert'
    for k=1:size(d,1)
      d{k,2} = ~d{k,2};
    end
end
set(h,'Data',d);
end
%===============================================================================

%===============================================================================
function showbaseline(src,~)
hFig = findall(0,'Tag','esfitFigure');
h = findobj(hFig,'Tag','baselinedata');
if src.Value
  set(h,'Visible','on')
else
  set(h,'Visible','off')
end
end
%===============================================================================

%===============================================================================
function showresiduals(src,~)
hFig = findall(0,'Tag','esfitFigure');
h(1) = findobj(hFig,'Tag','residualdata');
h(2) = findobj(hFig,'Tag','residualzero');
if src.Value
  set(h(1),'Visible','on')
  set(h(2),'Visible','on')
else
  set(h(1),'Visible','off')
  set(h(2),'Visible','off')
end
updateaxislimits()
end
%===============================================================================

%===============================================================================
function updateaxislimits(~)

global esfitdata
hFig = findall(0,'Tag','esfitFigure');

expdata = esfitdata.data(:);
if isfield(esfitdata,'curr') && isfield(esfitdata.curr,'sim')
  currsim = real(esfitdata.curr.sim(:));
else
  currsim = zeros(size(expdata));
end
mask = esfitdata.Opts.mask;
if isfield(esfitdata,'best') && isfield(esfitdata.best,'fit')
  bestsim = real(esfitdata.best.fit(:));
else
  bestsim = zeros(size(expdata));
end
plottedData = [expdata(mask); bestsim; currsim];
if get(findobj(hFig,'Tag','ResidualCheckbox'),'Value')
  residualzero = get(findobj(hFig,'Tag','residualzero'),'Value');
  residualdata = (expdata-currsim)+residualzero;
  plottedData = [plottedData; residualdata];
end
maxy = max(plottedData);
miny = min(plottedData);
YLimits = [miny maxy] + [-1 1]*esfitdata.Opts.PlotStretchFactor*(maxy-miny);

hAx = findobj(hFig,'Tag','dataaxes');
set(hAx,'YLim',YLimits);
x = get(findobj(hFig,'Tag','expdata'),'XData');
set(hAx,'XLim',[min(x) max(x)]);

end
%===============================================================================


%===============================================================================
function tableCellEditCallback(~,callbackData)
global esfitdata
hFig = findall(0,'Tag','esfitFigure');

% Get handle of table and row/column index of edited table cell
hTable = callbackData.Source;
ridx = callbackData.Indices(1);
cidx = callbackData.Indices(2);

if cidx==1
  allParamsFixed = all(~cell2mat(hTable.Data(:,1)));
  if allParamsFixed
    set(findobj(hFig,'Tag','StartButton'),'Enable','off');
  else
    set(findobj(hFig,'Tag','StartButton'),'Enable','on');
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
function axesButtonDownFcn(~,~)
global esfitdata
hFig = findall(0,'Tag','esfitFigure');
hAx = findobj(hFig,'Tag','dataaxes');

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
set(hFig,'WindowButtonMotionFcn',@(hObject,eventdata) drawmaskedregion(tmpmask));
set(hFig,'WindowButtonDownFcn',@(src,event)uiresume(src));
uiwait(hFig);
% waitforbuttonpress;
set(hFig,'WindowButtonMotionFcn',[])
set(hFig,'WindowButtonDownFcn',[])

% Update masked regions
cp = hAx.CurrentPoint;
x2 = cp(1,1);
maskrange = sort([x1 x2]);
esfitdata.Opts.mask(x>maskrange(1) & x<maskrange(2)) = 0;
delete(tmpmask);
showmaskedregions();
end
%===============================================================================

%===============================================================================
function drawmaskedregion(tmpmask)
hFig = findall(0,'Tag','esfitFigure');
hAx = findobj(hFig,'Tag','dataaxes');
cp = get(hAx,'CurrentPoint');
xdata = tmpmask.XData;
xdata(2:3) = cp(1,1);
set(tmpmask,'XData',xdata);
end
%===============================================================================

%===============================================================================
function showmaskedregions()
global esfitdata
hFig = findall(0,'Tag','esfitFigure');
hAx = findobj(hFig,'Tag','dataaxes');

% Delete existing mask patches
hMaskPatches = findobj(hAx,'Tag','maskPatch');
delete(hMaskPatches);

% Show masked-out regions
maskColor = [1 1 1]*0.95;
edges = find(diff([1; esfitdata.Opts.mask(:); 1]));
excludedRegions = reshape(edges,2,[]).';
upperlimit = numel(esfitdata.Opts.x);
excludedRegions(excludedRegions>upperlimit) = upperlimit;
excludedRegions = esfitdata.Opts.x(excludedRegions);

% Add a patch for each masked region
nMaskPatches = size(excludedRegions,1);
for r = 1:nMaskPatches
  x_patch = excludedRegions(r,[1 2 2 1]);
  y_patch = hAx.YLim([1 1 2 2]);
  patch(hAx,x_patch,y_patch,maskColor,'Tag','maskPatch','EdgeColor','none');
end

% Reorder so that mask patches are in the back
c = hAx.Children([nMaskPatches+1:end, 1:nMaskPatches]);
hAx.Children = c;
drawnow limitrate

end
%===============================================================================

%===============================================================================
function clearMaskCallback(~,~)
global esfitdata
esfitdata.Opts.mask = true(size(esfitdata.Opts.mask));
showmaskedregions();
esfitdata.best = [];
esfitdata.rmsdhistory = [];
esfitdata.besthistory.rmsd = [];
esfitdata.besthistory.par = [];

% Readjust vertical range
mask = esfitdata.Opts.mask;
expdata = esfitdata.data(:);
maxy = max(expdata(mask));
miny = min(expdata(mask));
YLimits = [miny maxy] + [-1 1]*esfitdata.Opts.PlotStretchFactor*(maxy-miny);
hFig = findall(0,'Tag','esfitFigure');
set(findobj(hFig,'Tag','dataaxes'),'YLim',YLimits);
drawnow

end
%===============================================================================

%===============================================================================
function updateLogBox(msg)
hFig = findall(0,'Tag','esfitFigure');

hLogBox = findobj(hFig,'Tag','LogBox');
txt = get(hLogBox,'Value');
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
set(hLogBox,'Value',txt)
drawnow limitrate
scroll(hLogBox,'bottom')

end
%===============================================================================

%===============================================================================
function copyLog(~,~)
% Copy log to clipboard
hFig = findall(0,'Tag','esfitFigure');
txt = get(findobj(hFig,'Tag','LogBox'),'Value');
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
global esfitdata
hFig = findall(0,'Tag','esfitFigure');

% Reactivate UI components
set(findobj(hFig,'Tag','SaveButton'),'Enable','off');

if isfield(esfitdata,'FitSets') && numel(esfitdata.FitSets)>0
  set(findobj(hFig,'Tag','selectStartPointButtonSelected'),'Enable','on');
  set(findobj(hFig,'Tag','exportSetButton'),'Enable','on');
  set(findobj(hFig,'Tag','deleteSetButton'),'Enable','on');
end

% Hide stop button, show start button
set(findobj(hFig,'Tag','StartButton'),'Text','Start fitting');

% Re-enable other buttons
set(findobj(hFig,'Tag','EvaluateButton'),'Enable','on');
set(findobj(hFig,'Tag','ResetButton'),'Enable','on');

% Re-enable listboxes
set(findobj(hFig,'Tag','AlgorithmMenu'),'Enable','on');
set(findobj(hFig,'Tag','TargetMenu'),'Enable','on');
set(findobj(hFig,'Tag','AutoScaleMenu'),'Enable','on');
set(findobj(hFig,'Tag','BaseLineMenu'),'Enable','on');

% Re-enable parameter table and its selection controls
set(findobj(hFig,'Tag','selectAllButton'),'Enable','on');
set(findobj(hFig,'Tag','selectNoneButton'),'Enable','on');
set(findobj(hFig,'Tag','selectInvButton'),'Enable','on');
set(findobj(hFig,'Tag','selectStartPointButtonCenter'),'Enable','on');
set(findobj(hFig,'Tag','selectStartPointButtonRandom'),'Enable','on');
set(findobj(hFig,'Tag','selectStartPointButtonBest'),'Enable','on');
colEditable = get(findobj(hFig,'Tag','ParameterTable'),'UserData');
set(findobj(hFig,'Tag','ParameterTable'),'ColumnEditable',colEditable);
set(findobj(hFig,'Tag','ParameterTable'),'CellEditCallback',@tableCellEditCallback);

% Re-enable mask tools
set(findobj(hFig,'Tag','clearMaskButton'),'Enable','on');
set(findobj(hFig,'Tag','MaskCheckbox'),'Enable','on');

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

global esfitdata
hPopup = findall(0,'Tag','algorithmpopup');

if strcmp(hPopup.Visible,'off')
  tabs = findobj(hPopup,'Tag','AlgorithmTabs');
  AlgorithmID = esfitdata.Opts.AlgorithmID;
  if AlgorithmID<=numel(tabs.Children)
    tabs.SelectedTab = tabs.Children(AlgorithmID);
  end
  set(findobj(hPopup,'Tag','nParamsField'),'Value',sum(~esfitdata.fixedParams));
  set(hPopup,'Visible','on')
elseif strcmp(hPopup.Visible,'on')
  set(hPopup,'Visible','off');
end

end
%===============================================================================

%===============================================================================
function selectAlgorithm(src)
% Callback for selection of fitting algorithm

global esfitdata

hFig = findall(0,'Tag','esfitFigure');
h(1) = findobj(hFig,'Tag','AlgorithmMenu');

hPopup = findall(0,'Tag','algorithmpopup');
h(2) = findobj(hPopup,'Tag','AlgorithmMenuPopup');
tabs = findobj(hPopup,'Tag','AlgorithmTabs');

% Update selection in main GUI window and algorithm settings popup
if nargin==1
  AlgorithmID = src.Value;
  h(1).Value = AlgorithmID;
  if AlgorithmID<=numel(tabs.Children)
    h(2).Value = AlgorithmID;
  end

  if AlgorithmID<=numel(tabs.Children)
    tabs.SelectedTab = tabs.Children(AlgorithmID);
  end
else
  AlgorithmID = h(1).Value;
end
esfitdata.Opts.AlgorithmID = AlgorithmID;

% Update the FitOpt structure
N = esfitdata.nParameters; %#ok<NASGU> 

if AlgorithmID<=numel(tabs.Children)
  hsetting = tabs.UserData{AlgorithmID};
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

global esfitdata
hPopup = findall(0,'Tag','algorithmpopup');
tabs = findobj(hPopup,'Tag','AlgorithmTabs');
nTabs = numel(tabs.UserData);

nParams = esfitdata.nParameters;

for i = 1:nTabs
  hsetting = tabs.UserData{i};
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

global esfitdata
N = esfitdata.nParameters; %#ok<NASGU> % needed to eval() expressions

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
hFig = findall(0,'Tag','esfitFigure');
AlgorithmID = get(findobj(hFig,'Tag','AlgorithmMenu'),'Value');

hPopup = findall(0,'Tag','algorithmpopup');
tabs = findobj(hPopup,'Tag','AlgorithmTabs');
IDselected = find(tabs.Children==tabs.SelectedTab);

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
function closeGUI(~,~)
% Callback to close GUI

global esfitdata;
esfitdata.UserCommand = 99;
drawnow;
hPopup = findall(0,'Tag','algorithmpopup');
delete(hPopup);
hFig = findall(0,'Tag','esfitFigure');
delete(hFig);

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

