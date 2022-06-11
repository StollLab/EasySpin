% esfit   Least-squares fitting for EPR and other data
%
%   esfit(data,fcn,p0,vary)
%   esfit(data,fcn,p0,lb,ub)
%   esfit(___,FitOpt)
%
%   pfit = esfit(___)
%   [pfit,datafit] = esfit(___)
%   [pfit,datafit,residuals] = esfit(___)
%
% Input:
%     data        experimental data, a vector of data points
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
%        .AutoScale either 1 (on) or 0 (off); default 1
%        .OutArg  two numbers [nOut iOut], where nOut is the number of
%                 outputs of the simulation function and iOut is the index
%                 of the output argument to use for fitting
%        .mask    array of 1 and 0 the same size as data vector
%                 values with mask 0 are excluded from the fit
% Output:
%     fit           structure with fitting results
%       .pfit       fitted parameter vector
%       .argsfit    fitted input arguments (if EasySpin-style)
%       .fitraw     fit, as returned by the simulation/model function
%       .fit        fit, including the fitted scale factor
%       .scale      fitted scale factor
%       .ci95       95% confidence intervals for all parameters
%       .corr       correlation matrix for all parameters
%       .cov        covariance matrix for all parameters
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
if ~isnumeric(data) || ~isvector(data) || isempty(data)
  error('First argument must be numeric experimental data in the form of a vector.');
end
esfitdata.data = data;

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
EasySpinFunction = any(strcmp(esfitdata.fcnName,{'pepper','garlic','chili','salt','curry'}));


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
esfitdata.pvec_lb = pvec_lb;
esfitdata.pvec_ub = pvec_ub;
esfitdata.nParameters = nParameters;

% Initialize parameter fixing mask used in GUI table
esfitdata.fixedParams = false(1,numel(pvec_0));

% Experimental parameters (for EasySpin functions)
%-------------------------------------------------------------------------------
if EasySpinFunction
  % Check or set Exp.nPoints
  if isfield(p0{2},'nPoints')
    if p0{2}.nPoints~=numel(data)
      error('Exp.nPoints is %d, but the data vector has %d elements.',...
        p0{2}.nPoints,numel(data));
    end
  else
    p0{2}.nPoints = numel(data);
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

keywords = strread(Opt.Method,'%s');
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
  Opt.mask = true(size(data));
else
  Opt.mask = logical(Opt.mask);
  if numel(Opt.mask)~=numel(data)
    error('Opt.mask has %d elements, but the data has %d elements.',numel(Opt.mask),numel(data));
  end
end
Opt.useMask = true;

% Scale fitting
if ~isfield(Opt,'AutoScale')
  Opt.AutoScale = 1;
end
switch Opt.AutoScale
  case 0, AutoScale = 0;
  case 1, AutoScale = 1;
  otherwise, error('Unknown setting for Opt.AutoScale - possible values are 0 and 1.');
end
esfitdata.AutoScale = AutoScale;
esfitdata.AutoScaleSettings = {1, 0};
esfitdata.AutoScaleStrings = {'on', 'off'};

StartpointNames{1} = 'center of range';
StartpointNames{2} = 'random within range';
StartpointNames{3} = 'selected parameter set';
esfitdata.StartpointNames = StartpointNames;

if ~isfield(Opt,'PrintLevel'), Opt.PrintLevel = 1; end

% Algorithm parameters
if ~isfield(Opt,'nTrials'), Opt.nTrials = 20000; end
if ~isfield(Opt,'TolFun'), Opt.TolFun = 1e-4; end
if ~isfield(Opt,'TolStep'), Opt.TolStep = 1e-6; end
if ~isfield(Opt,'maxTime'), Opt.maxTime = inf; end
if isfield(Opt,'RandomStart') && Opt.RandomStart
    Opt.Startpoint = 2; % random start point
else
    Opt.Startpoint = 1; % start point at center of range
end

% Grid search parameters
if ~isfield(Opt,'GridSize'), Opt.GridSize = 7; end
if ~isfield(Opt,'maxGridPoints'), Opt.maxGridPoints = 1e5; end
if ~isfield(Opt,'randomizeGrid'), Opt.randomizeGrid = true; end

% x axis for plotting
if ~isfield(Opt,'x')
  Opt.x = 1:numel(esfitdata.data);
end

% Internal parameters
if ~isfield(Opt,'PlotStretchFactor'), Opt.PlotStretchFactor = 0.05; end
if ~isfield(Opt,'maxParameters'), Opt.maxParameters = 30; end

if esfitdata.nParameters>Opt.maxParameters
    error('Cannot fit more than %d parameters simultaneously.',...
        Opt.maxParameters);
end
Opt.IterationPrintFunction = @iterationprint;

esfitdata.Opts = Opt;

esfitdata.maskSelectMode = false;

% Setup GUI and return if in interactive mode
%-------------------------------------------------------------------------------
interactiveMode = nargout==0;
if interactiveMode
  setupGUI(data);
  return
end

% Close GUI window if running in script mode
%-------------------------------------------------------------------------------
hFig = findobj('Tag','esfitFigure');
if ~isempty(hFig)
  close(hFig);
end

% Report parsed inputs
%-------------------------------------------------------------------------------
if esfitdata.Opts.PrintLevel>=1
  siz = size(esfitdata.data);
  if esfitdata.Opts.AutoScale
    autoScaleStr = 'on';
  else
    autoScaleStr = 'off';
  end
  fprintf('-- esfit ------------------------------------------------\n');
  fprintf('Data size:                [%d, %d]\n',siz(1),siz(2));
  fprintf('Model function name:      %s\n',esfitdata.fcnName);
  fprintf('Number of fit parameters: %d\n',esfitdata.nParameters);
  fprintf('Minimization algorithm:   %s\n',esfitdata.AlgorithmNames{esfitdata.Opts.AlgorithmID});
  fprintf('Residuals computed from:  %s\n',esfitdata.TargetNames{esfitdata.Opts.TargetID});
  fprintf('Autoscaling:              %s\n',autoScaleStr);
  fprintf('---------------------------------------------------------\n');
end

% Run least-squares fitting
%-------------------------------------------------------------------------------
result = runFitting();

% Report fit results
%-------------------------------------------------------------------------------
if esfitdata.Opts.PrintLevel>=1
  disp('---------------------------------------------------------');
  fprintf('Goodness of fit:\n');
  fprintf('   ssr             %g\n',result.ssr);
  fprintf('   rmsd            %g\n',result.rmsd);
  fprintf('   noise std       %g (estimated from residuals; assumes excellent fit)\n',std(result.residuals));
  fprintf('   chi-squared     %g (using noise std estimate; upper limit)\n',result.rmsd^2/var(result.residuals));
  if Opt.AutoScale
    fprintf('Fitted scale:       %g\n',result.scale);
  end
  fprintf('Parameters:\n');
  printparlist(result.pfit,esfitdata.pinfo,result.pstd,result.ci95);
  if ~isempty(result.corr) && numel(result.pfit)>1
    fprintf('Correlation matrix:\n');
    Sigma = result.corr;
    disp(Sigma);
    triuCorr = triu(abs(Sigma),1);
    fprintf('Strongest correlations:\n');
    [~,idx] = sort(triuCorr(:),'descend');
    [i1,i2] = ind2sub(size(Sigma),idx);
    np = numel(result.pfit);
    for k = 1:min(5,(np-1)*np/2)
      fprintf('    p(%d)-p(%d):    %g\n',i1(k),i2(k),Sigma(i1(k),i2(k)));
    end
    if any(reshape(triuCorr,1,[])>0.8)
      disp('    WARNING! Strong correlations between parameters.');
    end
  end
  disp('=========================================================');
end

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
nParameters = numel(esfitdata.pvec_0);
data_ = esfitdata.data;
fixedParams = esfitdata.fixedParams;
activeParams = ~fixedParams;
printLevel = esfitdata.Opts.PrintLevel;

% Set starting point
%-------------------------------------------------------------------------------
lb = esfitdata.pvec_lb;
ub = esfitdata.pvec_ub;
p0 = esfitdata.pvec_0;
switch esfitdata.Opts.Startpoint
  case 1 % provided start value
    p_start = p0;
  case 2 % random
    p_start = lb + rand(nParameters,1).*(ub-lb);
    p_start(fixedParams) = p0(fixedParams);
  case 3 % selected parameter set
    if useGUI
      h = findobj('Tag','SetListBox');
      s = h.String;
      if ~isempty(s)
        s = s{h.Value};
        ID = sscanf(s,'%d');
        idx = find([esfitdata.FitSets.ID]==ID);
        if ~isempty(idx)
          p_start = esfitdata.FitSets(idx).bestx;
        else
          error('Could not locate selected parameter set.');
        end
      else
        error('No saved parameter set yet.');
      end
    else
      error('Setting can only be used in GUI mode.');
    end
end
esfitdata.p_start = p_start;

% Run minimization over space of active parameters
%-------------------------------------------------------------------------------
fitOpt = esfitdata.Opts;
nActiveParams = sum(activeParams);
if nActiveParams>0
  if printLevel>=1
    fprintf('Running optimization algorithm with %d active parameters...\n',nActiveParams);
  end
  if useGUI
    fitOpt.IterFcn = @iterupdateGUI;
  end
  residualfun = @(x) residuals_(x,data_,fitOpt);
  rmsdfun = @(x) rmsd_(x,data_,fitOpt);
  p0_active = p_start(activeParams);
  lb_active = lb(activeParams);
  ub_active = ub(activeParams);
  switch fitOpt.AlgorithmID
    case 1 % Nelder-Mead simplex
      pfit_active = esfit_simplex(rmsdfun,p0_active,lb_active,ub_active,fitOpt);
    case 2 % Levenberg-Marquardt
      fitOpt.Gradient = fitOpt.TolFun;
      pfit_active = esfit_levmar(residualfun,p0_active,lb_active,ub_active,fitOpt);
    case 3 % Monte Carlo
      pfit_active = esfit_montecarlo(rmsdfun,lb_active,ub_active,fitOpt);
    case 4 % Genetic
      pfit_active = esfit_genetic(rmsdfun,lb_active,ub_active,fitOpt);
    case 5 % Grid search
      pfit_active = esfit_grid(rmsdfun,lb_active,ub_active,fitOpt);
    case 6 % Particle swarm
      pfit_active = esfit_swarm(rmsdfun,lb_active,ub_active,fitOpt);
    case 7 % lsqnonlin from Optimization Toolbox
      pfit_active = lsqnonlin(residualfun,p0_active,lb_active,ub_active);
  end
  pfit = p_start;
  pfit(activeParams) = pfit_active;
else
  if printLevel>=1
    disp('No active parameters; skipping optimization');
  end
  pfit = p_start;
end

if esfitdata.structureInputs
  argsfit = esfitdata.p2args(pfit);
else
  argsfit = [];
end
fit = esfitdata.bestfit;  % bestfit is set in residuals_
scale = esfitdata.bestscale;  % bestscale is set in residuals_
fitraw = fit/scale;

% Calculate metrics for goodness of fit
%-------------------------------------------------------------------------------
residuals = calculateResiduals(fit(:),esfitdata.data(:),esfitdata.Opts.TargetID);
residuals = residuals.'; % col -> row
ssr = sum(abs(residuals).^2); % sum of squared residuals
rmsd = sqrt(mean(residuals.^2)); % root-mean-square deviation

% Calculate parameter uncertainties
%-------------------------------------------------------------------------------
calculateUncertainties = esfitdata.UserCommand==0 && nActiveParams>0;
if calculateUncertainties
  if printLevel>=1
    disp('Calculating parameter uncertainties...');
    disp('  Estimating Jacobian...');
  end
  %maxRelStep = min((ub-pfit),(pfit-lb))./pfit;
  residualfun = @(x)residuals_(x,data_,fitOpt);
  J = jacobianest(residualfun,pfit_active);
  if ~any(isnan(J(:)))
    if printLevel>=1
      disp('  Calculating parameter covariance matrix...');
    end

    % Calculate covariance matrix and standard deviations
    covmatrix = hccm(J,residuals,'HC1');
    pstd = sqrt(diag(covmatrix));

    % Calculate confidence intervals
    norm_icdf = @(p)-sqrt(2)*erfcinv(2*p); % inverse of standard normal cdf
    ci = @(pctl)norm_icdf(1/2+pctl/2)*sqrt(diag(covmatrix));
    pctl = 0.95;
    ci95 = pfit_active + ci(pctl)*[-1 1];

    % Calculate correlation matrix
    if printLevel>=1
      disp('  Calculating parameter correlation matrix...');
    end
    Q = diag(diag(covmatrix).^(-1/2));
    corrmatrix = Q*covmatrix*Q;
  else
    if printLevel>=1
      disp('  NaN elements in Jacobian, cannot calculate parameter uncertainties.');
    end
    pstd = [];
    ci95 = [];
    covmatrix = [];
    corrmatrix = [];
  end
else
  if printLevel>=1
    disp('Fitting stopped by user. Skipping uncertainty quantification.');
  end
  pstd = [];
  ci95 = [];
  covmatrix = [];
  corrmatrix = [];
end

% Assemble output structure
%-------------------------------------------------------------------------------
result.argsfit = argsfit;
result.fit = fit;
result.scale = scale;
result.fitraw = fitraw;

result.pfit = pfit_active;
result.pnames = {esfitdata.pinfo.Name}.';
result.pnames = result.pnames(activeParams);

result.residuals = residuals;
result.ssr = ssr;
result.rmsd = rmsd;
result.pstd = pstd;
result.ci95 = ci95;
result.cov = covmatrix;
result.corr = corrmatrix;

end



%===============================================================================
function rmsd = rmsd_(x,data,Opt)
[~,rmsd] = residuals_(x,data,Opt);
end
%===============================================================================


%===============================================================================
function [residuals,rmsd,simdata,simscale] = residuals_(x,expdata,Opt)

% reads:
%   esfitdata.p_start
%   esfitdata.fixedParams
%   esfitdata.structureInputs
%   esfitdata.p2args
%   esfitdata.nOutArguments
%   esfitdata.OutArgument
%   esfitdata.fcn
% writes:
%   esfitdata.smallestError
%   esfitdata.bestfit
%   esfitdata.bestpar
%   esfitdata.errorlist

global esfitdata

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
if esfitdata.structureInputs
  args = esfitdata.p2args(par);
  try
    [out{:}] = esfitdata.fcn(args{:});
  catch ME
    error('\nThe model simulation function raised the following error:\n  %s\n',ME.message);
  end
else
  try
    [out{:}] = esfitdata.fcn(par);
  catch ME
    error('\nThe model simulation function raised the following error:\n  %s\n',ME.message);
  end
end
simdata = out{esfitdata.OutArgument}; % pick appropriate output argument

% Rescale simulated data if scale should be ignored
if Opt.AutoScale
  [~,simscale] = rescaledata(simdata(mask),expdata(mask),'lsq');
  simdata = simdata*simscale;
else
  simdata = simdata(:);
  simscale = 1;
end
esfitdata.currsim = simdata;
esfitdata.currpar = par;
esfitdata.currscale = simscale;

% Compute residuals
%-------------------------------------------------------------------------------
residuals = calculateResiduals(simdata(:),expdata(:),Opt.TargetID,mask(:));
rmsd = sqrt(mean(abs(residuals).^2));

% Keep track of errors
%-------------------------------------------------------------------------------
if ~isfield(esfitdata,'smallestError') || isempty(esfitdata.smallestError)
  esfitdata.smallestError = inf;
end
if ~isfield(esfitdata,'errorlist')
  esfitdata.errorlist = [];
end

esfitdata.errorlist = [esfitdata.errorlist rmsd];

isNewBest = rmsd<esfitdata.smallestError;
if isNewBest
  esfitdata.smallestError = rmsd;
  esfitdata.bestfit = simdata;
  esfitdata.bestscale = simscale;
  esfitdata.bestpar = par;
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
bestsim = real(esfitdata.bestfit(:));
currsim = real(esfitdata.currsim(:));
currpar = esfitdata.currpar;
bestx = currpar;
bestx(~esfitdata.fixedParams) = info.bestx;

% Update plotted data
set(findobj('Tag','expdata'),'XData',x,'YData',expdata);
set(findobj('Tag','bestsimdata'),'XData',x,'YData',bestsim);
set(findobj('Tag','currsimdata'),'XData',x,'YData',currsim);

% Readjust vertical range
mask = esfitdata.Opts.mask;
plottedData = [expdata(mask); bestsim; currsim];
maxy = max(plottedData);
miny = min(plottedData);
YLimits = [miny maxy] + [-1 1]*esfitdata.Opts.PlotStretchFactor*(maxy-miny);
set(findobj('Tag','dataaxes'),'YLim',YLimits);
drawnow

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
  oldvaluestring = striphtml(data{p,4});
  newvaluestring = sprintf('%0.6f',currpar(p));
  % Find first character at which the new value differs from the old one
  idx = 1;
  while idx<=min(length(oldvaluestring),length(newvaluestring))
    if oldvaluestring(idx)~=newvaluestring(idx), break; end
    idx = idx + 1;
  end
  active = data{p,1};
  if active
    str = ['<html><font color="#000000">' newvaluestring(1:idx-1) '</font><font color="#ff0000">' newvaluestring(idx:end) '</font></html>'];
  else
    str = ['<html><font color="#888888">' newvaluestring '</font></html>'];
  end
  data{p,4} = str;
end

% Update column with best values if current parameter set is new best
if info.newbest

  str = sprintf(' best RMSD: %g\n',esfitdata.smallestError);
  hRmsText = findobj('Tag','RmsText');
  set(hRmsText,'String',str);

  for p = 1:nParams
    oldvaluestring = striphtml(data{p,3});
    newvaluestring = sprintf('%0.6g',bestx(p));
    % Find first character at which the new value differs from the old one
    idx = 1;
    while idx<=min(length(oldvaluestring),length(newvaluestring))
      if oldvaluestring(idx)~=newvaluestring(idx), break; end
      idx = idx + 1;
    end
    active = data{p,1};
    if active
      str = ['<html><font color="#000000">' newvaluestring(1:idx-1) '</font><font color="#009900">' newvaluestring(idx:end) '</font></html>'];
    else
      str = ['<html><font color="#888888">' newvaluestring '</font></html>'];
    end
    data{p,3} = str;
  end
end

hParamTable.Data = data;

% Update rmsd plot
hErrorLine = findobj('Tag','errorline');
if ~isempty(hErrorLine)
  n = min(100,numel(esfitdata.errorlist));
  set(hErrorLine,'XData',1:n,'YData',log10(esfitdata.errorlist(end-n+1:end)));
  ax = hErrorLine.Parent;
  axis(ax,'tight');
  drawnow
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
  str = [indent sprintf('name%svalue        standard deviation        95%% confidence interval',repmat(' ',1,max(maxNameLength-4,0)+2))];
  for p = 1:nParams
    pname = pad(pinfo(p).Name,maxNameLength);
    str_ = sprintf('%s  %-#12.7g %-#12.7g (%6.3f %%)   %-#12.7g - %-#12.7g',pname,par(p),pstd(p),pstd(p)/par(p)*100,pci95(p,1),pci95(p,2));
    str = [str newline indent str_];
  end
else
  str = [indent sprintf('name%svalue',repmat(' ',1,max(maxNameLength-4,0)+2))];
  for p = 1:nParams
    pname = pad(pinfo(p).Name,maxNameLength);
    str_ = sprintf('%s  %-#12.7g',pname,par(p));
    str = [str newline indent str_];
  end
end

if nargout==0
  disp(str);
end

end
%===============================================================================


%===============================================================================
function residuals = calculateResiduals(A,B,mode,includemask)
AB = A - B;
if nargin>3
  AB(~includemask) = 0;
end
switch mode
  case 1  % fcn
    residuals = AB;
  case 2  % int
    residuals = cumsum(AB);
  case 3  % iint
    residuals = cumsum(cumsum(AB));
  case 4  % fft
    residuals = abs(fft(AB));
  case 5  % diff
    residuals = deriv(AB);
end
idxNaN = isnan(A) | isnan(B);
residuals(idxNaN) = 0; % ignore residual if A or B is NaN
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

% Disable listboxes
set(findobj('Tag','AlgorithMenu'),'Enable','off');
set(findobj('Tag','TargetMenu'),'Enable','off');
set(findobj('Tag','AutoScaleMenu'),'Enable','off');
set(findobj('Tag','StartpointMenu'),'Enable','off');

% Disable parameter table
set(findobj('Tag','selectAllButton'),'Enable','off');
set(findobj('Tag','selectNoneButton'),'Enable','off');
set(findobj('Tag','selectInvButton'),'Enable','off');
set(findobj('Tag','ParameterTable'),'Enable','off');

% Disable fitset list controls
set(findobj('Tag','deleteSetButton'),'Enable','off');
set(findobj('Tag','exportSetButton'),'Enable','off');
set(findobj('Tag','sortIDSetButton'),'Enable','off');
set(findobj('Tag','sortRMSDSetButton'),'Enable','off');

% Disable mask tools
set(findobj('Tag','clearMaskButton'),'Enable','off');
set(findobj('Tag','MaskMenu'),'Enable','off');

% Change GUI status
h = findobj('Tag','statusText');
if ~strcmp(h.String,'running')
  set(h,'String','running');
  drawnow
end

% Pull settings from UI
%-------------------------------------------------------------------------------
% Determine selected method, target, autoscaling, start point
esfitdata.Opts.AlgorithmID = get(findobj('Tag','AlgorithMenu'),'Value');
esfitdata.Opts.TargetID = get(findobj('Tag','TargetMenu'),'Value');
esfitdata.Opts.AutoScale = esfitdata.AutoScaleSettings{get(findobj('Tag','AutoScaleMenu'),'Value')};
esfitdata.Opts.Startpoint = get(findobj('Tag','StartpointMenu'),'Value');
esfitdata.Opts.useMask = get(findobj('Tag','MaskMenu'),'Value')==1;

% Get fixed parameters
data = get(findobj('Tag','ParameterTable'),'Data');
for iPar = 1:esfitdata.nParameters
  esfitdata.fixedParams(iPar) = data{iPar,1}==0;
end

% Run fitting
%-------------------------------------------------------------------------------
useGUI = true;
result = runFitting(useGUI);

% Save result to fit set list
esfitdata.currFitSet = result;
esfitdata.currFitSet.Target = esfitdata.TargetNames{esfitdata.Opts.TargetID};


% Update GUI with fit results
%-------------------------------------------------------------------------------

set(findobj('Tag','statusText'),'String','');

% Remove current values from parameter table
hTable = findobj('Tag','ParameterTable');
Data = hTable.Data;
for p = 1:size(Data,1), Data{p,4} = '-'; end
set(hTable,'Data',Data);

% Hide current sim plot in data axes
set(findobj('Tag','currsimdata'),'YData',NaN(1,numel(esfitdata.data)));
hErrorLine = findobj('Tag','errorline');
set(hErrorLine,'XData',1,'YData',NaN);
axis(hErrorLine.Parent,'tight');
drawnow
set(findobj('Tag','logLine'),'String','');

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

% Re-enable listboxes
set(findobj('Tag','AlgorithMenu'),'Enable','on');
set(findobj('Tag','TargetMenu'),'Enable','on');
set(findobj('Tag','AutoScaleMenu'),'Enable','on');
set(findobj('Tag','StartpointMenu'),'Enable','on');

% Re-enable parameter table and its selection controls
set(findobj('Tag','selectAllButton'),'Enable','on');
set(findobj('Tag','selectNoneButton'),'Enable','on');
set(findobj('Tag','selectInvButton'),'Enable','on');
set(findobj('Tag','ParameterTable'),'Enable','on');

% Re-enable mask tools
set(findobj('Tag','clearMaskButton'),'Enable','on');
set(findobj('Tag','MaskMenu'),'Enable','on');

end
%===============================================================================


%===============================================================================
function iterationprint(str)
hLogLine = findobj('Tag','logLine');
if isempty(hLogLine)
  disp(str);
else
  set(hLogLine,'String',str);
end
end
%===============================================================================


%===============================================================================
function str = striphtml(str)
html = false;
for k = numel(str):-1:1
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
    values = fitset.pfit;

    % Set column with best-fit parameter values
    hTable = findobj('Tag','ParameterTable');
    data = get(hTable,'data');
    for p = 1:numel(values)
      data{p,3} = sprintf('%0.6g',values(p));
    end
    set(hTable,'Data',data);

    h = findobj('Tag','bestsimdata');
    set(h,'YData',fitset.fit);
    drawnow
  end
else
  h = findobj('Tag','bestsimdata');
  set(h,'YData',NaN(size(h.YData)));
  drawnow
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
fitSet = rmfield(esfitdata.FitSets(idx),'bestx');
varname = sprintf('fit%d',ID);
assignin('base',varname,fitSet);
fprintf('Fit set %d assigned to variable ''%s''.\n',ID,varname);
evalin('base',varname);
end
%===============================================================================


%===============================================================================
function selectAllButtonCallback(~,~)
h = findobj('Tag','ParameterTable');
d = h.Data;
d(:,1) = {true};
set(h,'Data',d);
end
%===============================================================================


%===============================================================================
function selectNoneButtonCallback(~,~)
h = findobj('Tag','ParameterTable');
d = h.Data;
d(:,1) = {false};
set(h,'Data',d);
end
%===============================================================================


%===============================================================================
function selectInvButtonCallback(~,~)
h = findobj('Tag','ParameterTable');
d = h.Data;
for k=1:size(d,1)
    d{k,1} = ~d{k,1};
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
  s{k} = sprintf('%d. rmsd %g (%s)',...
    esfitdata.FitSets(k).ID,esfitdata.FitSets(k).rmsd,esfitdata.FitSets(k).Target);
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
function clearMaskCallback(~,~)
global esfitdata
esfitdata.Opts.mask = true(size(esfitdata.Opts.mask));
showmaskedregions();
end
%===============================================================================


%===============================================================================
function saveFitsetCallback(~,~)
global esfitdata
esfitdata.lastSetID = esfitdata.lastSetID+1;
esfitdata.currFitSet.ID = esfitdata.lastSetID;
if ~isfield(esfitdata,'FitSets') || isempty(esfitdata.FitSets)
  esfitdata.FitSets(1) = esfitdata.currFitSet;
else
  esfitdata.FitSets(end+1) = esfitdata.currFitSet;
end
refreshFitsetList(-1);
end
%===============================================================================


%===============================================================================
function tableCellEditCallback(~,callbackData)
global esfitdata

% Get handle of table and row/column index of edited table cell
hTable = callbackData.Source;
ridx = callbackData.Indices(1);
cidx = callbackData.Indices(2);

% Return unless it's a cell that contains lower or upper bound
lbColumn = 5; % lower-bound column
ubColumn = 6; % upper-bound column
lbedit = cidx==lbColumn;
ubedit = cidx==ubColumn;
if ~lbedit && ~ubedit, return; end

% Get lower and upper bounds of interval from table
lower = str2double(hTable.Data{ridx,lbColumn});
upper = str2double(hTable.Data{ridx,ubColumn});

% Convert user-entered string to number
newval = str2double(callbackData.EditData);

% Revert if conversion didn't yield a scalar
if numel(newval)~=1 || isnan(newval) || ~isreal(newval)
  warning('Input ''%s'' is not a number. Reverting edit.',callbackData.EditData);
  hTable.Data{ridx,cidx} = callbackData.PreviousData;
  return
end

% Set new lower/upper bound
if lbedit
  lower = newval;
elseif ubedit
  upper = newval;
end

% Revert if lower bound would be above upper bound
if lower>upper
  warning('Lower bound is above upper bound. Reverting edit.');
  hTable.Data{ridx,cidx} = callbackData.PreviousData;
  return
end

% Update lower and upper bounds
esfitdata.pvec_lb(ridx) = lower;
esfitdata.pvec_ub(ridx) = upper;

end
%===============================================================================


%===============================================================================
function setupGUI(data)

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

sz = [1400 800]; % figure size
screensize = get(0,'ScreenSize');
xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the screen vertically
set(hFig,'position',[xpos, ypos, sz(1), sz(2)],'units','pixels');
set(hFig,'WindowStyle','normal','DockControls','off','MenuBar','none');
set(hFig,'Resize','off');
set(hFig,'Name','esfit - Least-Squares Fitting','NumberTitle','off');
set(hFig,'CloseRequestFcn',...
    'global esfitdata; esfitdata.UserCommand = 99; drawnow; delete(gcf);');

% Axes
%-------------------------------------------------------------------------------
% data display
hAx = axes('Parent',hFig,'Tag','dataaxes','Units','pixels',...
    'Position',[30 30 1010 740],'FontSize',8,'Layer','top');
x0 = 1060; % Start of display to the right of the axes

NaNdata = NaN(1,numel(data));
mask = esfitdata.Opts.mask;
dispData = esfitdata.data;
maxy = max(dispData(mask));
miny = min(dispData(mask));
YLimits = [miny maxy] + [-1 1]*Opt.PlotStretchFactor*(maxy-miny);
minx = min(esfitdata.Opts.x);
maxx = max(esfitdata.Opts.x);
x = esfitdata.Opts.x;

h(1) = line(x,NaNdata,'Color','k','Marker','.','LineStyle','none');
h(2) = line(x,NaNdata,'Color',[0 0.6 0]);
h(3) = line(x,NaNdata,'Color','r');
set(h(1),'Tag','expdata','XData',esfitdata.Opts.x,'YData',dispData);
set(h(2),'Tag','bestsimdata');
set(h(3),'Tag','currsimdata');
hAx.XLim = [minx maxx];
hAx.YLim = YLimits;
hAx.ButtonDownFcn = @axesButtonDownFcn;
grid(hAx,'on');
%set(hAx,'XTick',[],'YTick',[]);
box on

showmaskedregions();

% iteration and rms error displays
%-------------------------------------------------------------------------------
y0 = 160;
hAx = axes('Parent',hFig,'Units','pixels','Position',[x0 y0 100 80],'Layer','top');
h = plot(hAx,1,NaN,'.');
set(h,'Tag','errorline','MarkerSize',5,'Color',[0.2 0.2 0.8]);
set(gca,'FontSize',7,'YScale','lin','XTick',[],'YAxisLoc','right','Layer','top','YGrid','on');
title('log10(RMSD)','Color','k','FontSize',7,'FontWeight','normal');

h = uicontrol('Style','text','Position',[x0+125 y0+64 205 16]);
set(h,'FontSize',8,'String',' RMSD: -','ForegroundColor',[0 0 1],'Tooltip','Current best RMSD');
set(h,'Tag','RmsText','HorizontalAl','left');

h = uicontrol('Style','text','Position',[x0+125 y0 205 62]);
set(h,'FontSize',7,'Tag','logLine','Tooltip','Information from fitting algorithm');
set(h,'Horizontal','left');


% Parameter table
%-------------------------------------------------------------------------------
columnname = {'','Name','best','current','lower','upper'};
columnformat = {'logical','char','char','char','char','char'};
colEditable = [true false false false true true];
data = cell(numel(esfitdata.pinfo),6);
for p = 1:numel(esfitdata.pinfo)
  data{p,1} = true;
  data{p,2} = char(esfitdata.pinfo(p).Name);
  data{p,3} = '-';
  data{p,4} = '-';
  data{p,5} = sprintf('%0.6g',esfitdata.pvec_lb(p));
  data{p,6} = sprintf('%0.6g',esfitdata.pvec_ub(p));
end
y0 = 500; dx = 80;
uitable('Tag','ParameterTable',...
    'FontSize',8,...
    'Position',[x0 y0 330 250],...
    'ColumnFormat',columnformat,...
    'ColumnName',columnname,...
    'ColumnEditable',colEditable,...
    'CellEditCallback',@tableCellEditCallback,...
    'ColumnWidth',{20,62,62,62,62,60},...
    'RowName',[],...
    'Data',data);
uicontrol('Style','text',...
    'Position',[x0 y0+250 230 20],...
    'BackgroundColor',get(gcf,'Color'),...
    'FontWeight','bold','String','Parameters',...
    'HorizontalAl','left');
uicontrol('Style','pushbutton','Tag','selectInvButton',...
    'Position',[x0+210 y0+250 50 20],...
    'String','invert','Enable','on','Callback',@selectInvButtonCallback,...
    'HorizontalAl','left',...
    'Tooltip','Invert selection of parameters');
uicontrol('Style','pushbutton','Tag','selectAllButton',...
    'Position',[x0+260 y0+250 30 20],...
    'String','all','Enable','on','Callback',@selectAllButtonCallback,...
    'HorizontalAl','left',...
    'Tooltip','Select all parameters');
uicontrol('Style','pushbutton','Tag','selectNoneButton',...
    'Position',[x0+290 y0+250 40 20],...
    'String','none','Enable','on','Callback',@selectNoneButtonCallback,...
    'HorizontalAl','left',...
    'Tooltip','Unselect all parameters');
uicontrol(hFig,'Style','text',...
    'String','Function',...
    'Tooltip','Name of simulation function',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0+270 dx 20]);
uicontrol(hFig,'Style','text','Tag','statusText',...
    'String','',...
    'Tooltip','status',...
    'FontWeight','bold',...
    'ForegroundColor','r',...
    'HorizontalAlign','right',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0+270 y0+270 60 20]);
uicontrol(hFig,'Style','text',...
    'String',esfitdata.fcnName,...
    'ForeGroundColor','b',...
    'Tooltip',sprintf('using output no. %d of %d',esfitdata.nOutArguments,esfitdata.OutArgument),...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0+dx y0+270 dx 20]);

% popup menus
%-------------------------------------------------------------------------------
dx = 60; y0 = 390; dy = 24;
uicontrol(hFig,'Style','text',...
    'String','Algorithm',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0+3*dy-4 dx 20]);
uicontrol(hFig,'Style','popupmenu',...
    'Tag','AlgorithMenu',...
    'String',esfitdata.AlgorithmNames,...
    'Value',Opt.AlgorithmID,...
    'BackgroundColor','w',...
    'Tooltip','Fitting algorithm',...
    'Position',[x0+dx y0+3*dy 150 20]);

uicontrol(hFig,'Style','text',...
    'String','Target',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0+2*dy-4 dx 20]);
uicontrol(hFig,'Style','popupmenu',...
    'Tag','TargetMenu',...
    'String',esfitdata.TargetNames,...
    'Value',Opt.TargetID,...
    'BackgroundColor','w',...
    'Tooltip','Target function',...
    'Position',[x0+dx y0+2*dy 150 20]);

uicontrol(hFig,'Style','text',...
    'String','AutoScale',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0+dy-4 dx 20]);
uicontrol(hFig,'Style','popupmenu',...
    'Tag','AutoScaleMenu',...
    'String',esfitdata.AutoScaleStrings,...
    'Value',find(cellfun(@(x)x==esfitdata.AutoScale,esfitdata.AutoScaleSettings),1),...
    'BackgroundColor','w',...
    'Tooltip','Autoscaling',...
    'Position',[x0+dx y0+dy 150 20]);

uicontrol(hFig,'Style','text',...
    'String','Startpoint',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0-4 dx 20]);
h = uicontrol(hFig,'Style','popupmenu',...
    'Tag','StartpointMenu',...
    'String',esfitdata.StartpointNames,...
    'Value',1,...
    'BackgroundColor','w',...
    'Tooltip','Starting point for fit',...
    'Position',[x0+dx y0 150 20]);
if esfitdata.Opts.Startpoint==2, set(h,'Value',2); end

uicontrol(hFig,'Style','text',...
    'String','Use mask',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0-4-dy dx 20]);
uicontrol(hFig,'Style','popupmenu',...
    'Tag','MaskMenu',...
    'String',{'on','off'},...
    'Value',1,...
    'BackgroundColor','w',...
    'Tooltip','Use mask',...
    'Position',[x0+dx y0-dy 150 20]);


% Start/Stop buttons
%--------------------------------------------------------------------------
pos =  [x0+220 y0-3+30 110 67];
pos1 = [x0+220 y0-3    110 30];
pos2 = [x0+220 y0-3-30    110 30];
uicontrol(hFig,'Style','pushbutton',...
    'Tag','StartButton',...
    'String','Start',...
    'Callback',@startButtonCallback,...
    'Visible','on',...
    'Tooltip','Start fitting',...
    'Position',pos);
uicontrol(hFig,'Style','pushbutton',...
    'Tag','StopButton',...
    'String','Stop',...
    'Visible','off',...
    'Tooltip','Stop fitting',...
    'Callback','global esfitdata; esfitdata.UserCommand = 1;',...
    'Position',pos);
uicontrol(hFig,'Style','pushbutton',...
    'Tag','SaveButton',...
    'String','Save parameter set',...
    'Callback',@saveFitsetCallback,...
    'Enable','off',...
    'Tooltip','Save latest fitting result',...
    'Position',pos1);
uicontrol(hFig,'Style','pushbutton',...
    'Tag','SaveButton',...
    'String','Clear mask',...
    'Callback',@clearMaskCallback,...
    'Enable','on',...
    'Tooltip','Clear mask',...
    'Position',pos2);

% Fitset list
%-------------------------------------------------------------------------------
y0 = 10;
uicontrol('Style','text','Tag','SetListTitle',...
    'Position',[x0 y0+100 230 20],...
    'BackgroundColor',get(gcf,'Color'),...
    'FontWeight','bold','String','Parameter sets',...
    'Tooltip','List of stored fit parameter sets',...
    'HorizontalAl','left');
uicontrol(hFig,'Style','listbox','Tag','SetListBox',...
    'Position',[x0 y0 330 100],...
    'String','','Tooltip','',...
    'BackgroundColor',[1 1 0.9],...
    'KeyPressFcn',@deleteSetListKeyPressFcn,...
    'Callback',@setListCallback);
uicontrol(hFig,'Style','pushbutton','Tag','deleteSetButton',...
    'Position',[x0+280 y0+100 50 20],...
    'String','delete',...
    'Tooltip','Delete fit set','Enable','off',...
    'Callback',@deleteSetButtonCallback);
uicontrol(hFig,'Style','pushbutton','Tag','exportSetButton',...
    'Position',[x0+230 y0+100 50 20],...
    'String','export',...
    'Tooltip','Export fit set to workspace','Enable','off',...
    'Callback',@exportSetButtonCallback);
uicontrol(hFig,'Style','pushbutton','Tag','sortIDSetButton',...
    'Position',[x0+210 y0+100 20 20],...
    'String','id',...
    'Tooltip','Sort parameter sets by ID','Enable','off',...
    'Callback',@sortIDSetButtonCallback);
uicontrol(hFig,'Style','pushbutton','Tag','sortRMSDSetButton',...
    'Position',[x0+180 y0+100 30 20],...
    'String','rmsd',...
    'Tooltip','Sort parameter sets by rmsd','Enable','off',...
    'Callback',@sortRMSDSetButtonCallback);

drawnow

set(hFig,'NextPlot','new');

end

function axesButtonDownFcn(~,~)
global esfitdata
hAx = findobj('Tag','dataaxes');
cp = hAx.CurrentPoint;
x = esfitdata.Opts.x;
maskSelectMode = esfitdata.maskSelectMode;
if maskSelectMode
  x1 = esfitdata.maskSelect.x1;
  x2 = cp(1,1);
  maskrange = sort([x1 x2]);
  esfitdata.Opts.mask(x>maskrange(1) & x<maskrange(2)) = 0;
  showmaskedregions();
else
  esfitdata.maskSelect.x1 = cp(1,1);
end
esfitdata.maskSelectMode = ~esfitdata.maskSelectMode;
end


function showmaskedregions()
global esfitdata
hAx = findobj('Tag','dataaxes');

% Delete existing mask patches
hMaskPatches = findobj(hAx,'Tag','maskPatch');
delete(hMaskPatches);

% Show masked-out regions
maskColor = [1 1 1]*0.95;
edges = find(diff([1; esfitdata.Opts.mask(:); 1]));
excludedRegions = reshape(edges,2,[]).';
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
drawnow

end
