% esfit   Least-squares fitting for EPR data
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
%     p0          starting input parameters
%                   EasySpin-style functions: {Sys0,Exp0} or {Sys0,Exp0,Opt0}
%                   other functions: vector
%     vary        allowed variation of parameters
%                   EasySpin-style functions: {vSys} or {vSys,vExp} or {vSys,vExp,vOpt}
%                   other functions: vector
%     lb          lower bounds of parameters
%                   EasySpin-style functions: {lbSys,lbExp} or {lbSys,lbExp,lbOpt}
%                   other functions: vector
%     ub          upper bounds of parameters
%                   EasySpin-style functions: {ubSys,ubExp} or {ubSys,ubExp,ubOpt}
%                   other functions: vector
%     FitOpt      options for the fitting algorithms
%        .Method  string containing kewords for
%           -algorithm: 'simplex','levmar','montecarlo','genetic','grid','swarm'
%           -target function: 'fcn', 'int', 'dint', 'diff', 'fft'
%        .Scaling string with scaling method keyword
%           'maxabs' (default), 'lsq', 'lsq0','lsq1','lsq2','none'
%        .OutArg  two numbers [nOut iOut], where nOut is the number of
%                 outputs of the simulation function and iOut is the index
%                 of the output argument to use for fitting
% Output:
%     fit           structure with fitting results
%       .pfit       fitted parameter vector
%       .argsfit    fitter input arguments (if EasySpin-style)
%       .sim        simulated data
%       .simscaled  simulated data, scaled to the experimental ones
%       .ci95       95% confidence intervals for all parameters
%       .corr       correlation matrix for all parameters
%       .cov        covariance matrix for all parameters
%

function result = esfit(data,fcn,p0,varargin)

if nargin==0
    
    Sys.g = 2.00123;
    Sys.lwpp = [0.0 0.8];
    Exp.mwFreq = 9.5;
    Exp.Range = [330 350];
    [B,spc] = pepper(Sys,Exp);
    rng(123415);
    Ampl = 100;
    spc = Ampl*addnoise(spc,50,'n');
    
    Sys0.g = 2.003;
    Sys0.lwpp = [0.0 0.7];
    vSys.g = 0.01;
    vSys.lwpp = [0.0 0.25];
    
    Opt = struct;
    FitOpt = struct;
    FitOpt.Method = 'levmar fcn';
    FitOpt.TolFun = 1e-6;
    FitOpt.PrintLevel = 2;
    fit = esfit(spc,@pepper,{Sys0,Exp,Opt},{vSys},FitOpt)
    
    subplot(2,1,1)
    plot(B,spc,B,fit.sim);
    subplot(2,1,2)
    plot(B,spc-fit.sim);
    return
end

if nargin==0, help(mfilename); return; end

% Check expiry date
error(eschecker);

if nargin<4
  error('At least 3 inputs are required (data, fcn, p0, pvary).');
end
if nargin>6
  error('At most 6 inputs are accepted.');
end

% Parse argument list
varyProvided = true;
switch nargin
  case 4
    pvary = varargin{1};
    FitOpt = struct;
  case 5
    if isstruct(varargin{2})
      pvary = varargin{1};
      FitOpt = varargin{2};
    else
      varyProvided = false;
      lb = varargin{1};
      ub = varargin{2};
      FitOpt = struct;
    end
  case 6
    varyProvided = false;
    lb = varargin{1};
    ub = varargin{2};
    FitOpt = varargin{3};
end

if isempty(FitOpt)
  FitOpt = struct;
end
if ~isstruct(FitOpt)
  error('FitOpt (last input argument) must be a structure.');
end

% Set up global structure for data sharing among local functions
global fitdat
fitdat.currFitSet = [];

% Load utility functions
argspar = esfit_argsparams();


% Experimental data
%-------------------------------------------------------------------------------
if isstruct(data) || ~isnumeric(data)
    error('First argument must be numeric experimental data.');
end
fitdat.nSpectra = 1;
fitdat.ExpSpec = data;
fitdat.ExpSpecScaled = rescale(data,'maxabs');


% Model function
%-------------------------------------------------------------------------------
if ~isa(fcn,'function_handle')
  str = 'The simulation/model function (2nd input) must be a function handle (with @).';
  if ischar(fcn)
    error('%s\nUse esfit(data,@%s,...) instead of esfit(data,''%s'',...).',str,fcn,fcn);
  else
    error('%s\nFor example, to use the function pepper(...), use esfit(data,@pepper,...).',str);
  end
end
try
  nargin(fcn);
catch
  error('The simulation/model function given as second input cannot be found.');
end

fitdat.fcn = fcn;
fitdat.fcnName = func2str(fcn);

fitdat.lastSetID = 0;

EasySpinFunction = any(strcmp(fitdat.fcnName,{'pepper','garlic','chili','salt','curry'}));


% Parameters
%-------------------------------------------------------------------------------
structureInputs = iscell(p0);
fitdat.structureInputs = structureInputs;

% Determine parameter intervals, either from p0 and pvary, or from lower/upper bounds
if structureInputs
  argspar.validargs(p0);
  fitdat.nSystems = numel(p0{1});
  if varyProvided
    % use p0 and pvary
    pinfo = argspar.getparaminfo(pvary);
    argspar.checkparcompatibility(pinfo,p0);
    pvec_0 = argspar.getparamvalues(p0,pinfo);
    pvec_vary = argspar.getparamvalues(pvary,pinfo);
    pvec_lb = pvec_0 - pvec_vary;
    pvec_ub = pvec_0 + pvec_vary;
  else
    % use lower and upper bounds
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
  for k = numel(p0):-1:1
    pinfo(k).Name = sprintf('p(%d)',k);
  end
end

% Assure all parameter vectors are column vectors
pvec_0 = pvec_0(:);
pvec_lb = pvec_lb(:);
pvec_ub = pvec_ub(:);

fitdat.args = p0;
fitdat.pinfo = pinfo;
fitdat.pvec_0 = pvec_0;
fitdat.pvec_lb = pvec_lb;
fitdat.pvec_ub = pvec_ub;
fitdat.nParameters = numel(pvec_0);
if fitdat.nParameters==0
    error('No variable parameters to fit.');
end
fitdat.inactiveParams = false(numel(pvec_0),1);
fitdat.p2args = @(pars) argspar.setparamvalues(p0,pinfo,pars);


% Experimental parameters (for EasySpin functions)
%-------------------------------------------------------------------------------
if EasySpinFunction
  % Set Exp.nPoints
  if isfield(p0{2},'nPoints')
    if p0{2}.nPoints~=numel(data)
      error('Exp.nPoints is %d, but the spectral data vector is %d long.',...
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


% Fitting options
%===============================================================================
if ~isfield(FitOpt,'OutArg')
  fitdat.nOutArguments = abs(nargout(fitdat.fcn));
  fitdat.OutArgument = fitdat.nOutArguments;
else
  if numel(FitOpt.OutArg)~=2
    error('FitOpt.OutArg must contain two values [nOut iOut]');
  end
  if FitOpt.OutArg(2)>FitOpt.OutArg(1)
    error('FitOpt.OutArg: second number cannot be larger than first one.');
  end
  fitdat.nOutArguments = FitOpt.OutArg(1);
  fitdat.OutArgument = FitOpt.OutArg(2);  
end

if ~isfield(FitOpt,'Scaling'), FitOpt.Scaling = 'lsq0'; end

if ~isfield(FitOpt,'Method'), FitOpt.Method = ''; end
FitOpt.MethodID = 1; % simplex
FitOpt.TargetID = 1; % function as is
if EasySpinFunction
  if isfield(p0{2},'Harmonic') && p0{2}.Harmonic>0
    FitOpt.TargetID = 2; % integral
  else
    if strcmp(fitdat.fcnName,'pepper') || strcmp(fitdat.fcnName,'garlic')
      FitOpt.TargetID = 2; % integral
    end
  end
end

keywords = strread(FitOpt.Method,'%s'); %#ok<DSTRRD>
for k = 1:numel(keywords)
  switch keywords{k}
    case 'simplex',    FitOpt.MethodID = 1;
    case 'levmar',     FitOpt.MethodID = 2;
    case 'montecarlo', FitOpt.MethodID = 3;
    case 'genetic',    FitOpt.MethodID = 4;
    case 'grid',       FitOpt.MethodID = 5;
    case 'swarm',      FitOpt.MethodID = 6;
    case 'lsqnonlin',  FitOpt.MethodID = 7;
      
    case 'fcn',        FitOpt.TargetID = 1;
    case 'int',        FitOpt.TargetID = 2;
    case 'iint',       FitOpt.TargetID = 3;
    case 'dint',       FitOpt.TargetID = 3;
    case 'diff',       FitOpt.TargetID = 4;
    case 'fft',        FitOpt.TargetID = 5;
    otherwise
      error('Unknown ''%s'' in FitOpt.Method.',keywords{k});
  end
end

MethodNames{1} = 'Nelder-Mead simplex';
MethodNames{2} = 'Levenberg-Marquardt';
MethodNames{3} = 'Monte Carlo';
MethodNames{4} = 'genetic algorithm';
MethodNames{5} = 'grid search';
MethodNames{6} = 'particle swarm';
MethodNames{7} = 'lsqnonlin';
fitdat.MethodNames = MethodNames;

TargetNames{1} = 'data as is';
TargetNames{2} = 'integral';
TargetNames{3} = 'double integral';
TargetNames{4} = 'derivative';
TargetNames{5} = 'Fourier transform';
fitdat.TargetNames = TargetNames;

ScalingNames{1} = 'scale only (max abs)';
ScalingNames{2} = 'scale only (lsq)';
ScalingNames{3} = 'scale & shift (lsq0)';
ScalingNames{4} = 'scale & linear baseline (lsq1)';
ScalingNames{5} = 'scale & quad. baseline (lsq2)';
ScalingNames{6} = 'no scaling';
fitdat.ScalingNames = ScalingNames;

ScalingString{1} = 'maxabs';
ScalingString{2} = 'lsq';
ScalingString{3} = 'lsq0';
ScalingString{4} = 'lsq1';
ScalingString{5} = 'lsq2';
ScalingString{6} = 'none';
fitdat.ScalingString = ScalingString;

StartpointNames{1} = 'center of range';
StartpointNames{2} = 'random within range';
StartpointNames{3} = 'selected parameter set';
fitdat.StartpointNames = StartpointNames;

FitOpt.ScalingID = find(strcmp(FitOpt.Scaling,ScalingString));
if isempty(FitOpt.ScalingID)
    error('Unknown ''%s'' in FitOpt.Scaling.',FitOpt.Scaling);
end

fitdat.GUI = nargout==0;

if ~isfield(FitOpt,'PrintLevel'), FitOpt.PrintLevel = 1; end
if ~isfield(FitOpt,'nTrials'), FitOpt.nTrials = 20000; end
if ~isfield(FitOpt,'TolFun'), FitOpt.TolFun = 1e-4; end
if ~isfield(FitOpt,'TolStep'), FitOpt.TolStep = 1e-6; end
if ~isfield(FitOpt,'maxTime'), FitOpt.maxTime = inf; end
if isfield(FitOpt,'RandomStart') && FitOpt.RandomStart
    FitOpt.Startpoint = 2; % random start point
else
    FitOpt.Startpoint = 1; % start point at center of range
end

if ~isfield(FitOpt,'GridSize'), FitOpt.GridSize = 7; end

% Internal parameters
if ~isfield(FitOpt,'PlotStretchFactor'), FitOpt.PlotStretchFactor = 0.05; end
if ~isfield(FitOpt,'maxGridPoints'), FitOpt.maxGridPoints = 1e5; end
if ~isfield(FitOpt,'maxParameters'), FitOpt.maxParameters = 30; end
if fitdat.nParameters>FitOpt.maxParameters
    error('Cannot fit more than %d parameters simultaneously.',...
        FitOpt.maxParameters);
end
FitOpt.IterationPrintFunction = @iterationprint;
fitdat.FitOpts = FitOpt;

% Setup UI and quit if in GUI mode
if fitdat.GUI
  setupGUI(data);
  return
end

% Run least-squares fitting if not in GUI mode
%[fit_args,fit_sim,fit_residuals,~,fit_par] = runFitting();
out = runFitting();

% Collect result output structure
if structureInputs
  result.afit = out.argsfit;
end
result.pfit = out.pfit;
result.rmsd = out.rmsd;
result.ci95 = out.ci95;
result.corr = out.corr;
result.cov = out.cov;
result.residuals = out.residuals;
result.sim = out.fitSpec;
result.simscaled = out.fitSpecScaled;

clear global UserCommand

end

%===============================================================================
%===============================================================================
%===============================================================================


%===============================================================================
% Run fitting algorithm
%===============================================================================
function result = runFitting()

global UserCommand fitdat
UserCommand = 0;

if fitdat.FitOpts.PrintLevel
  fprintf('-- esfit ------------------------------------------------\n');
  fprintf('Function name:            %s\n',fitdat.fcnName);
  fprintf('Number of parameters:     %d\n',fitdat.nParameters);
  fprintf('Number of datasets:       %d\n',fitdat.nSpectra);
  fprintf('Minimization method:      %s\n',fitdat.MethodNames{fitdat.FitOpts.MethodID});
  fprintf('Residuals computed from:  %s\n',fitdat.TargetNames{fitdat.FitOpts.TargetID});
  fprintf('Scaling mode:             %s\n',fitdat.FitOpts.Scaling);
  fprintf('---------------------------------------------------------\n');
end

% Set starting point
%-------------------------------------------------------------------------------
lb = fitdat.pvec_lb;
ub = fitdat.pvec_ub;
p_center = (lb+ub)/2;
switch fitdat.FitOpts.Startpoint
  case 1 % center of range
    p_start = p_center;
  case 2 % random
    p_start = lb + rand(fitdat.nParameters,1).*(ub-lb);
    p_start(fitdat.inactiveParams) = p_center(fitdat.inactiveParams);
  case 3 % selected parameter set
    h = findobj('Tag','SetListBox');
    s = h.String;
    if ~isempty(s)
      s = s{h.Value};
      ID = sscanf(s,'%d');
      idx = find([fitdat.FitSets.ID]==ID);
      if ~isempty(idx)
        p_start = fitdat.FitSets(idx).bestx;
      else
        error('Could not locate selected parameter set.');
      end
    else
      error('No saved parameter set yet.');
    end
end
fitdat.p_start = p_start;

if strcmp(fitdat.FitOpts.Scaling,'none')
    data_ = fitdat.ExpSpec;
else
    data_ = fitdat.ExpSpecScaled;
end

% Run minimization over space of active parameters
%-------------------------------------------------------------------------------
pfit = p_start;
active = ~fitdat.inactiveParams;
fitOpts = fitdat.FitOpts;
if sum(active)>0
  residualfun = @(x)residuals_(x,data_,fitdat,fitOpts);
  rmsdfun = @(x)rmsd_(x,data_,fitdat,fitOpts);
  p0_a = p_start(active);
  lb_a = lb(active);
  ub_a = ub(active);
  switch fitOpts.MethodID
    case 1 % Nelder-Mead simplex
      pfit_a = esfit_simplex(rmsdfun,p0_a,lb_a,ub_a,fitOpts);
    case 2 % Levenberg-Marquardt
      fitOpts.Gradient = fitOpts.TolFun;
      pfit_a = esfit_levmar(residualfun,p0_a,lb_a,ub_a,fitOpts);
    case 3 % Monte Carlo
      pfit_a = esfit_montecarlo(rmsdfun,lb_a,ub_a,fitOpts);
    case 4 % Genetic
      pfit_a = esfit_genetic(rmsdfun,lb_a,ub_a,fitOpts);
    case 5 % Grid search
      pfit_a = esfit_grid(rmsdfun,lb_a,ub_a,fitOpts);
    case 6 % Particle swarm
      pfit_a = esfit_swarm(rmsdfun,lb_a,ub_a,fitOpts);
    case 7 % lsqnonlin from Optimization Toolbox
      pfit_a = lsqnonlin(residualfun,p0_a,lb_a,ub_a);
  end
  pfit(active) = pfit_a;
end

% Finish up
%-------------------------------------------------------------------------------

% Simulate model fit
if fitdat.structureInputs
  argsfit = fitdat.p2args(pfit);
  [out{1:fitdat.nOutArguments}] = fitdat.fcn(argsfit{:});
else
  argsfit = [];
  [out{1:fitdat.nOutArguments}] = fitdat.fcn(pfit);
end

fitSpec = out{fitdat.OutArgument}; % pick last output argument

fitSpecScaled = rescale(fitSpec,fitdat.ExpSpecScaled,fitdat.FitOpts.Scaling);
fitSpec =       rescale(fitSpec,fitdat.ExpSpec,      fitdat.FitOpts.Scaling);

% Calculate residuals and rmsd
residuals = calculateResiduals(fitSpecScaled(:),fitdat.ExpSpecScaled(:),fitdat.FitOpts.TargetID);
residuals = residuals.'; % col -> row
rmsd = sqrt(mean(residuals.^2));

% Calculate parameter confidence intervals and correlation matrix
%-------------------------------------------------------------------------------
disp('Calculating parameter uncertainties...');
disp('  Estimating Jacobian...');
maxRelStep = min((ub-pfit),(pfit-lb))./pfit;
J = jacobianest(residualfun,pfit,maxRelStep);
if any(isnan(J(:)))
  disp('  NaN elements in Jacobian, cannot calculate parameter uncertainties.');
  covmatrix = [];
  corrmatrix = [];
  ci = @(pctl)[];
  UQdone = false;
else
  disp('  Calculating parameter covariance matrix...');
  pinv(J.'*J)
  covmatrix = hccm(J,residuals,'HC1');
  
  % Calculate correlation matrix
  disp('  Calculating parameter correlation matrix...');
  Q = diag(diag(covmatrix).^(-1/2));
  corrmatrix = Q*covmatrix*Q;
  
  norm_icdf = @(p)-sqrt(2)*erfcinv(2*p); % inverse of standard normal cdf
  ci = @(pctl)norm_icdf(1/2+pctl/2)*sqrt(diag(covmatrix));
  UQdone = true;
end

% Report
%-------------------------------------------------------------------------------
if fitdat.FitOpts.PrintLevel && UserCommand~=99
  disp('---------------------------------------------------------');
  fprintf('Goodness of fit:\n');
  fprintf('   ssr             %g\n',sum(residuals.^2));
  fprintf('   rmsd            %g\n',rmsd);
  fprintf('   noise std       %g (estimated from residuals)\n',std(residuals));
  fprintf('   chi^2           %g (using noise std estimate; upper limit)\n',rmsd^2/var(residuals));
  if UQdone
    pctl = 0.95;
    ci95 = pfit + ci(pctl)*[-1 1];
    fprintf('Best-fit parameter values (%d%% confidence intervals):\n',100*pctl);
  else
    disp('Best-fit parameter values without confidence intervals:');
    ci95 = [];
  end
  str = printparlist(pfit,fitdat.pinfo,ci95);
  fprintf(str);
  if ~isempty(corrmatrix)
    fprintf('Correlation matrix:\n');
    disp(corrmatrix);
    triuCorr = triu(abs(corrmatrix),1);
    fprintf('Strongest correlations:\n');
    [~,idx] = sort(triuCorr(:),'descend');
    [i1,i2] = ind2sub(size(corrmatrix),idx);
    np = numel(pfit);
    for k = 1:min(5,(np-1)*np/2)
      fprintf('    p(%d)-p(%d):    %g\n',i1(k),i2(k),corrmatrix(i1(k),i2(k)));
    end
    if any(reshape(triuCorr,1,[])>0.8)
      disp('    WARNING! Stong correlations between parameters.');
    end
  end
  disp('=========================================================');
end

% Assemble output
%-------------------------------------------------------------------------------
result.argsfit = argsfit;
result.fitSpec = fitSpec;
result.residuals = residuals;
result.fitSpecScaled = fitSpecScaled;
result.pfit = pfit;
result.ci95 = ci95;
result.cov = covmatrix;
result.corr = corrmatrix;
result.rmsd = rmsd;

end

%===============================================================================
function startButton_cb(~,~)

global fitdat UserCommand

UserCommand = 0;

% Update UI, pull settings from UI
%-------------------------------------------------------------------------------

% Hide Start button, show Stop button
set(findobj('Tag','StopButton'),'Visible','on');
set(findobj('Tag','StartButton'),'Visible','off');
set(findobj('Tag','SaveButton'),'Enable','off');

% Disable listboxes
set(findobj('Tag','MethodMenu'),'Enable','off');
set(findobj('Tag','TargetMenu'),'Enable','off');
set(findobj('Tag','ScalingMenu'),'Enable','off');
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

% Determine selected method, target, and scaling
fitdat.FitOpts.MethodID = get(findobj('Tag','MethodMenu'),'Value');
fitdat.FitOpts.TargetID = get(findobj('Tag','TargetMenu'),'Value');
fitdat.FitOpts.Scaling = fitdat.ScalingString{get(findobj('Tag','ScalingMenu'),'Value')};
fitdat.FitOpts.Startpoint = get(findobj('Tag','StartpointMenu'),'Value');

% Get inactive parameters
data = get(findobj('Tag','ParameterTable'),'Data');
for iPar = 1:fitdat.nParameters
    fitdat.inactiveParams(iPar) = data{iPar,1}==0;
end

% Change GUI status
h = findobj('Tag','statusText');
if ~strcmp(h.String,'running')
  set(h,'String','running');
  drawnow
end

% Run fitting
%-------------------------------------------------------------------------------
out = runFitting();

argsfit = out.argsfit;
BestSpec = out.fitSpec;
BestSpecScaled = out.fitSpecScaled;
bestx = out.pfit;

out.argsfit = argsfit;
out.fitSpec = BestSpec;
out.fitSpecScaled = BestSpecScaled;
out.pfit = xfit;

rmsd = sqrt(mean(residuals.^2));

% GUI update
%-------------------------------------------------------------------------------

set(findobj('Tag','statusText'),'String','');

% Remove current values from parameter table
hTable = findobj('Tag','ParameterTable');
Data = hTable.Data;
for p = 1:size(Data,1), Data{p,4} = '-'; end
set(hTable,'Data',Data);

% Hide current sim plot in data axes
set(findobj('Tag','currsimdata'),'YData',NaN*ones(1,numel(fitdat.ExpSpec)));
hErrorLine = findobj('Tag','errorline');
set(hErrorLine,'XData',1,'YData',NaN);
axis(hErrorLine.Parent,'tight');
drawnow
set(findobj('Tag','logLine'),'String','');

% Reactivate UI components
set(findobj('Tag','SaveButton'),'Enable','on');

if isfield(fitdat,'FitSets') && numel(fitdat.FitSets)>0
    set(findobj('Tag','deleteSetButton'),'Enable','on');
    set(findobj('Tag','exportSetButton'),'Enable','on');
    set(findobj('Tag','sortIDSetButton'),'Enable','on');
    set(findobj('Tag','sortRMSDSetButton'),'Enable','on');
end

% Hide stop button, show start button
set(findobj('Tag','StopButton'),'Visible','off');
set(findobj('Tag','StartButton'),'Visible','on');

% Re-enable listboxes
set(findobj('Tag','MethodMenu'),'Enable','on');
set(findobj('Tag','TargetMenu'),'Enable','on');
set(findobj('Tag','ScalingMenu'),'Enable','on');
set(findobj('Tag','StartpointMenu'),'Enable','on');

% Re-enable parameter table and its selection controls
set(findobj('Tag','selectAllButton'),'Enable','on');
set(findobj('Tag','selectNoneButton'),'Enable','on');
set(findobj('Tag','selectInvButton'),'Enable','on');
set(findobj('Tag','ParameterTable'),'Enable','on');

% Save current set to set list
newFitSet.rmsd = rmsd;
if strcmp(fitdat.FitOpts.Scaling,'none')
    newFitSet.fitSpec = BestSpec;
    newFitSet.expSpec = fitdat.ExpSpec;
else
    newFitSet.fitSpec = BestSpecScaled;
    newFitSet.expSpec = fitdat.ExpSpecScaled;
end
newFitSet.residuals = residuals;
newFitSet.bestx = bestx;
%newFitSet.bestvalues = bestvalues;
TargetKey = {'fcn','int','iint','diff','fft'};
newFitSet.Target = TargetKey{fitdat.FitOpts.TargetID};
if numel(argsfit)==1
    newFitSet.Sys = argsfit{1};
else
    newFitSet.Sys = argsfit;
end
fitdat.currFitSet = newFitSet;

end
%===============================================================================


%===============================================================================
function rmsd = rmsd_(x,ExpSpec,FitData,FitOpt)
[~,rmsd] = residuals_(x,ExpSpec,FitData,FitOpt);
end
%===============================================================================


%===============================================================================
function varargout = residuals_(x,expdata,FitInfo,FitOpt)

global UserCommand

if ~isfield(FitInfo,'smallestError') || isempty(FitInfo.smallestError)
    FitInfo.smallestError = inf;
end
if ~isfield(FitInfo,'errorlist')
    FitInfo.errorlist = [];
end

% Evaluate model function-------------------------------------------------------
p_all = FitInfo.p_start;
active = ~FitInfo.inactiveParams;
p_all(active) = x;
par = p_all;
if FitInfo.structureInputs
  args = FitInfo.p2args(par);
  try
    [out{1:FitInfo.nOutArguments}] = FitInfo.fcn(args{:});
  catch ME
    error('\nThe model simulation function raised the following error:\n  %s\n',ME.message);
  end
else
  try
    [out{1:FitInfo.nOutArguments}] = FitInfo.fcn(par);
  catch ME
    error('\nThe model simulation function raised the following error:\n  %s\n',ME.message);
  end
end
simdata = out{FitInfo.OutArgument}; % pick appropriate output argument

% Scale simulated spectrum to experimental spectrum ----------------------------
simdata = rescale(simdata,expdata,FitOpt.Scaling);

% Compute residuals ------------------------------------------------------------
residuals = calculateResiduals(simdata(:),expdata(:),FitOpt.TargetID);
rmsd = real(sqrt(mean(residuals.^2)));

FitInfo.errorlist = [FitInfo.errorlist rmsd];
isNewBest = rmsd<FitInfo.smallestError;

if isNewBest
    FitInfo.smallestError = rmsd;
    FitInfo.bestspec = simdata;
    FitInfo.bestpar = par;
end

% Update GUI
%-------------------------------------------------------------------------------
if FitInfo.GUI && UserCommand~=99
    
    % update plot
    x = 1:numel(expdata);
    set(findobj('Tag','expdata'),'XData',x,'YData',expdata);
    set(findobj('Tag','bestsimdata'),'XData',x,'YData',real(FitInfo.bestspec));
    set(findobj('Tag','currsimdata'),'XData',x,'YData',real(simdata));
    
    % readjust vertical range
    dispData = [expdata(:); real(FitInfo.bestspec(:)); real(simdata(:))];
    maxy = max(dispData);
    miny = min(dispData);
    YLimits = [miny maxy] + [-1 1]*FitOpt.PlotStretchFactor*(maxy-miny);
    set(findobj('Tag','dataaxes'),'YLim',YLimits);
    drawnow
    
    % update numbers parameter table
    if UserCommand~=99
        
        % update column with current parameter values
        hParamTable = findobj('Tag','ParameterTable');
        data = get(hParamTable,'data');
        for p = 1:numel(par)
            olddata = striphtml(data{p,4});
            newdata = sprintf('%0.6f',par(p));
            idx = 1;
            while idx<=length(olddata) && idx<=length(newdata)
                if olddata(idx)~=newdata(idx), break; end
                idx = idx + 1;
            end
            active = data{p,1};
            if active
                data{p,4} = ['<html><font color="#000000">' newdata(1:idx-1) '</font><font color="#ff0000">' newdata(idx:end) '</font></html>'];
            else
                data{p,4} = ['<html><font color="#888888">' newdata '</font></html>'];
            end
        end
        
        % update column with best values if current parameter set is new best
        if isNewBest
            
            str = sprintf(' best RMSD: %g\n',(FitInfo.smallestError));
            hRmsText = findobj('Tag','RmsText');
            set(hRmsText,'String',str);
            
            for p = 1:numel(par)
                olddata = striphtml(data{p,3});
                newdata = sprintf('%0.6g',par(p));
                idx = 1;
                while idx<=length(olddata) && idx<=length(newdata)
                    if olddata(idx)~=newdata(idx), break; end
                    idx = idx + 1;
                end
                active = data{p,1};
                if active
                    data{p,3} = ['<html><font color="#000000">' newdata(1:idx-1) '</font><font color="#009900">' newdata(idx:end) '</font></html>'];
                else
                    data{p,3} = ['<html><font color="#888888">' newdata '</font></html>'];
                end
            end
        end
        set(hParamTable,'Data',data);
        
    end
    
    hErrorLine = findobj('Tag','errorline');
    if ~isempty(hErrorLine)
        n = min(100,numel(FitInfo.errorlist));
        set(hErrorLine,'XData',1:n,'YData',log10(FitInfo.errorlist(end-n+1:end)));
        ax = hErrorLine.Parent;
        axis(ax,'tight');
        drawnow
    end
    
end

if UserCommand==2
    UserCommand = 0;
    str = printparlist(pfit,fitdat.pinfo);
    disp('--- current best fit parameters -------------')
    fprintf(str);
    disp('---------------------------------------------')
end

out = {residuals,rmsd,simdata};
varargout = out(1:nargout);

end
%===============================================================================


%===============================================================================
% Print parameters.
function str = printparlist(par,pinfo,pci)
nParams = numel(par);
str = '';
maxNameLength = max(arrayfun(@(x)length(x.Name),pinfo));
if nargin==3 && ~isempty(pci)
  for p = 1:nParams
    str = sprintf('%s   %s  %g  (%g - %g)\n',str,pad(pinfo(p).Name,...
      maxNameLength),par(p),pci(p,1),pci(p,2));
  end
else
  for p = 1:nParams
    str = sprintf('%s   %s  %g\n',str,pad(pinfo(p).Name,maxNameLength),par(p));
  end
end

if nargout==0, fprintf(str); end
end
%===============================================================================


%===============================================================================
function residuals = calculateResiduals(A,B,mode)
switch mode
    case 1 % fcn
        residuals = A-B;
    case 2 % int
        residuals = cumsum(A-B);
    case 3 % iint
        residuals = cumsum(cumsum(A-B));
    case 4 % fft
        residuals = abs(fft(A-B));
    case 5 % diff
        residuals = deriv(A-B);
end
idxNaN = isnan(A) | isnan(B);
residuals(idxNaN) = 0; % ignore residual if A or B is NaN
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
global fitdat
h = findobj('Tag','SetListBox');
idx = h.Value;
str = h.String;
nSets = numel(str);
if nSets>0
    ID = sscanf(str{idx},'%d');
    for k = numel(fitdat.FitSets):-1:1
        if fitdat.FitSets(k).ID==ID
            fitdat.FitSets(k) = [];
        end
    end
    if idx>length(fitdat.FitSets), idx = length(fitdat.FitSets); end
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
global FitResults
h = findobj('Tag','SetListBox');
idx = h.Value;
str = h.String;
if ~isempty(str)
    ID = sscanf(str{idx},'%d');
    k = find([FitResults.FitSets.ID]==ID);
    if k>0
        fitset = FitResults.FitSets(k);
        values = fitset.bestx;
        
        % Set column with best-fit parameter values
        hTable = findobj('Tag','ParameterTable');
        data = get(hTable,'data');
        for p = 1:numel(values)
            data{p,3} = sprintf('%0.6g',values(p));
        end
        set(hTable,'Data',data);
        
        h = findobj('Tag','bestsimdata');
        set(h,'YData',fitset.fitSpec);
        drawnow
    end
else
    h = findobj('Tag','bestsimdata');
    set(h,'YData',h.YData*NaN);
    drawnow
end

end
%===============================================================================


%===============================================================================
function exportSetButtonCallback(~,~)
global FitResults
h = findobj('Tag','SetListBox');
v = h.Value;
s = h.String;
ID = sscanf(s{v},'%d');
idx = [FitResults.FitSets.ID]==ID;
fitSet = rmfield(FitResults.FitSets(idx),'bestx');
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
global FitResults
ID = [FitResults.FitSets.ID];
[~,idx] = sort(ID);
FitResults.FitSets = FitResults.FitSets(idx);
refreshFitsetList(0);
end
%===============================================================================


%===============================================================================
function sortRMSDSetButtonCallback(~,~)
global FitResults
rmsd = [FitResults.FitSets.rmsd];
[~,idx] = sort(rmsd);
FitResults.FitSets = FitResults.FitSets(idx);
refreshFitsetList(0);
end
%===============================================================================


%===============================================================================
function refreshFitsetList(idx)
global FitResults
h = findobj('Tag','SetListBox');
nSets = numel(FitResults.FitSets);
s = cell(1,nSets);
for k = 1:nSets
    s{k} = sprintf('%d. rmsd %g (%s)',...
        FitResults.FitSets(k).ID,FitResults.FitSets(k).rmsd,FitResults.FitSets(k).Target);
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
function saveFitsetCallback(~,~)
global FitResults
FitResults.lastSetID = FitResults.lastSetID+1;
FitResults.currFitSet.ID = FitResults.lastSetID;
if ~isfield(FitResults,'FitSets') || isempty(FitResults.FitSets)
    FitResults.FitSets(1) = FitResults.currFitSet;
else
    FitResults.FitSets(end+1) = FitResults.currFitSet;
end
refreshFitsetList(-1);
end
%===============================================================================


%===============================================================================
function tableCellEditCallback(~,callbackData)
global fitdat

hTable = callbackData.Source;

% Get row and column index of edited table cell
ridx = callbackData.Indices(1);
cidx = callbackData.Indices(2);

% Return unless it's the last two columns
if cidx<5, return; end

% Revert if user-entered string does not cleanly convert to a scalar.
newval = str2double(callbackData.EditData);
if numel(newval)~=1 || isnan(newval)
    hTable.Data{ridx,cidx} = callbackData.PreviousData;
    warning('Input is not a number.');
    return
end

% Get lower and upper bounds of interval from table
if cidx==5
    lower = newval;
    upper = str2double(hTable.Data{ridx,6});
elseif cidx==6
    lower = str2double(hTable.Data{ridx,5});
    upper = newval;
end

% Revert if lower bound would be above upper bound
if lower>upper
    warning('Lower bound is above upper bound.');
    hTable.Data{ridx,cidx} = callbackData.PreviousData;
    return
end

% Get parameter string (e.g. 'g(1)', or 'B.g(2)' for more than 1 system)
% and determine system index
parName = hTable.Data{ridx,2};
if fitdat.nSystems>1
    parName = parName(3:end); % disregard 'A.' etc.
    iSys = parName(1)-64; % 'A' -> 1, 'B' -> 2, etc
else
    iSys = 1;
end

% Update appropriate field in FitData.Sys0 and FitData.Vary
center = (lower+upper)/2;
vary = (upper-lower)/2;
eval(sprintf('FitData.Sys0{%d}.%s = %g;',iSys,parName,center));
eval(sprintf('FitData.Vary{%d}.%s = %g;',iSys,parName,vary));

end
%===============================================================================


%===============================================================================
function setupGUI(data)

global fitdat
FitOpt = fitdat.FitOpts;

% main figure
%---------------------------------------------------------------------------
hFig = findobj('Tag','esfitFigure');
if isempty(hFig)
    hFig = figure('Tag','esfitFigure','WindowStyle','normal');
else
    figure(hFig);
    clf(hFig);
end

sz = [1000 600]; % figure size
screensize = get(0,'ScreenSize');
xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the screen vertically
set(hFig,'position',[xpos, ypos, sz(1), sz(2)],'units','pixels');
set(hFig,'WindowStyle','normal','DockControls','off','MenuBar','none');
set(hFig,'Resize','off');
set(hFig,'Name','EasySpin Least-Squares Fitting','NumberTitle','off');
set(hFig,'CloseRequestFcn',...
    'global UserCommand; UserCommand = 99; drawnow; delete(gcf);');

% axes
%---------------------------------------------------------------------------
excludedRegions = [];
% data display
hAx = axes('Parent',hFig,'Units','pixels',...
    'Position',[10 10 640 580],'FontSize',8,'Layer','top');
NaNdata = ones(1,numel(data))*NaN;
dispData = fitdat.ExpSpecScaled;
maxy = max(dispData); miny = min(dispData);
YLimits = [miny maxy] + [-1 1]*FitOpt.PlotStretchFactor*(maxy-miny);
for r = 1:size(excludedRegions,1)
    h = patch(excludedRegions(r,[1 2 2 1]),YLimits([1 1 2 2]),[1 1 1]*0.8);
    set(h,'EdgeColor','none');
end
x = 1:numel(data);
h(1) = line(x,NaNdata,'Color','k','Marker','.');
h(2) = line(x,NaNdata,'Color',[0 0.6 0]);
h(3) = line(x,NaNdata,'Color','r');
set(h(1),'Tag','expdata','XData',1:numel(dispData),'YData',dispData);
set(h(2),'Tag','bestsimdata');
set(h(3),'Tag','currsimdata');
set(hAx,'XLim',[1 numel(dispData)]);
set(hAx,'YLim',YLimits);
set(hAx,'Tag', 'dataaxes');
set(hAx,'XTick',[],'YTick',[]);
box on

% iteration and rms error displays
%-----------------------------------------------------------------
x0 = 660; y0 = 160;
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
%---------------------------------------------------------------------------
columnname = {'','Name','best','current','lower','upper'};
columnformat = {'logical','char','char','char','char','char'};
colEditable = [true false false false true true];
data = cell(numel(fitdat.pinfo),6);
for p = 1:numel(fitdat.pinfo)
    data{p,1} = true;
    data{p,2} = char(fitdat.pinfo(p).Name);
    data{p,3} = '-';
    data{p,4} = '-';
    data{p,5} = sprintf('%0.6g',fitdat.pvec_lb(p));
    data{p,6} = sprintf('%0.6g',fitdat.pvec_ub(p));
end
x0 = 660; y0 = 400; dx = 80;
uitable('Tag','ParameterTable',...
    'FontSize',8,...
    'Position',[x0 y0 330 150],...
    'ColumnFormat',columnformat,...
    'ColumnName',columnname,...
    'ColumnEditable',colEditable,...
    'CellEditCallback',@tableCellEditCallback,...
    'ColumnWidth',{20,62,62,62,62,60},...
    'RowName',[],...
    'Data',data);
uicontrol('Style','text',...
    'Position',[x0 y0+150 230 20],...
    'BackgroundColor',get(gcf,'Color'),...
    'FontWeight','bold','String','Parameters',...
    'HorizontalAl','left');
uicontrol('Style','pushbutton','Tag','selectInvButton',...
    'Position',[x0+210 y0+150 50 20],...
    'String','invert','Enable','on','Callback',@selectInvButtonCallback,...
    'HorizontalAl','left',...
    'Tooltip','Invert selection of parameters');
uicontrol('Style','pushbutton','Tag','selectAllButton',...
    'Position',[x0+260 y0+150 30 20],...
    'String','all','Enable','on','Callback',@selectAllButtonCallback,...
    'HorizontalAl','left',...
    'Tooltip','Select all parameters');
uicontrol('Style','pushbutton','Tag','selectNoneButton',...
    'Position',[x0+290 y0+150 40 20],...
    'String','none','Enable','on','Callback',@selectNoneButtonCallback,...
    'HorizontalAl','left',...
    'Tooltip','Unselect all parameters');
uicontrol(hFig,'Style','text',...
    'String','Function',...
    'Tooltip','Name of simulation function',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0+170 dx 20]);
uicontrol(hFig,'Style','text','Tag','statusText',...
    'String','',...
    'Tooltip','status',...
    'FontWeight','bold',...
    'ForegroundColor','r',...
    'HorizontalAlign','right',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0+270 y0+170 60 20]);
uicontrol(hFig,'Style','text',...
    'String',fitdat.fcnName,...
    'ForeGroundColor','b',...
    'Tooltip',sprintf('using output no. %d of %d',fitdat.nOutArguments,fitdat.OutArgument),...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0+dx y0+170 dx 20]);

% popup menus
%---------------------------------------------------------------------------
x0 = 660; dx = 60; y0 = 290; dy = 24;
uicontrol(hFig,'Style','text',...
    'String','Method',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0+3*dy-4 dx 20]);
uicontrol(hFig,'Style','popupmenu',...
    'Tag','MethodMenu',...
    'String',fitdat.MethodNames,...
    'Value',FitOpt.MethodID,...
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
    'String',fitdat.TargetNames,...
    'Value',FitOpt.TargetID,...
    'BackgroundColor','w',...
    'Tooltip','Target function',...
    'Position',[x0+dx y0+2*dy 150 20]);
uicontrol(hFig,'Style','text',...
    'String','Scaling',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0+dy-4 dx 20]);
uicontrol(hFig,'Style','popupmenu',...
    'Tag','ScalingMenu',...
    'String',fitdat.ScalingNames,...
    'Value',FitOpt.ScalingID,...
    'BackgroundColor','w',...
    'Tooltip','Scaling mode',...
    'Position',[x0+dx y0+dy 150 20]);
uicontrol(hFig,'Style','text',...
    'String','Startpoint',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0-4 dx 20]);
h = uicontrol(hFig,'Style','popupmenu',...
    'Tag','StartpointMenu',...
    'String',fitdat.StartpointNames,...
    'Value',1,...
    'BackgroundColor','w',...
    'Tooltip','Starting point for fit',...
    'Position',[x0+dx y0 150 20]);
if fitdat.FitOpts.Startpoint==2, set(h,'Value',2); end

% Start/Stop buttons
%---------------------------------------------------------------------------
pos =  [x0+220 y0-3+30 110 67];
pos1 = [x0+220 y0-3    110 30];
uicontrol(hFig,'Style','pushbutton',...
    'Tag','StartButton',...
    'String','Start',...
    'Callback',@startButton_cb,...
    'Visible','on',...
    'Tooltip','Start fitting',...
    'Position',pos);
uicontrol(hFig,'Style','pushbutton',...
    'Tag','StopButton',...
    'String','Stop',...
    'Visible','off',...
    'Tooltip','Stop fitting',...
    'Callback','global UserCommand; UserCommand = 1;',...
    'Position',pos);
uicontrol(hFig,'Style','pushbutton',...
    'Tag','SaveButton',...
    'String','Save parameter set',...
    'Callback',@saveFitsetCallback,...
    'Enable','off',...
    'Tooltip','Save latest fitting result',...
    'Position',pos1);

% Fitset list
%---------------------------------------------------------------------------
x0 = 660; y0 = 10;
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
