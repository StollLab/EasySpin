% esfit   least-squares fitting for EPR spectral simulations
%
%   esfit(simfunc,expspc,Sys0,Vary,Exp)
%   esfit(simfunc,expspc,Sys0,Vary,Exp,SimOpt)
%   esfit(simfunc,expspc,Sys0,Vary,Exp,SimOpt,FitOpt)
%   bestsys = esfit(...)
%   [bestsys,bestspc] = esfit(...)
%
%     simfunc   simulation function name:
%                 'pepper', 'garlic', 'salt', 'chili', etc.
%     expspc    experimental spectrum, a vector of data points
%     Sys0      starting values for spin system parameters
%     Vary      allowed variation of parameters
%     Exp       experimental parameter, for simulation function
%     SimOpt    options for the simulation algorithms
%
%     FitOpt      options for the fitting algorithms
%        Method   string containing kewords for
%          -algorithm: 'simplex','levmar','montecarlo','genetic','grid'
%          -rms basis: 'fcn', 'int', 'dint', 'diff', 'fft'
%        Scaling  string with scaling method keyword
%          'maxabs' (default), 'minmax', 'lsq', 'lsq0','lsq1','lsq2'

function varargout = esfit(SimFunctionName,ExpSpec,Sys0,Vary,Exp,SimOpt,FitOpt)

if (nargin==0), help(mfilename); return; end

% --------License ------------------------------------------------
LicErr = 'Could not determine license.';
Link = 'epr@eth'; eschecker; error(LicErr); clear Link LicErr
% --------License ------------------------------------------------

if (nargin<5), error('Not enough inputs.'); end
if (nargin<6), SimOpt = struct('unused',NaN); end
if (nargin<7), FitOpt = struct('unused',NaN); end
  
global SimFcnName hAxes hText smallerr bestspec
hAxes = [];
hText = [];
smallerr = [];
%bestspec = [];

% Simulation function name
%--------------------------------------------------------------------
if ~ischar(SimFunctionName) && ~isa(SimFunctionName,'function_handle')
  error('First parameter must be simulation function name.');
end
SimFcnName = SimFunctionName;  

% System structure
%--------------------------------------------------------------------
if ~iscell(Sys0), Sys0 = {Sys0}; end
nSystems = numel(Sys0);
for s = 1:nSystems
  if ~isfield(Sys0{s},'weight'), Sys0{s}.weight = 1; end
end

% Experimental spectrum
%--------------------------------------------------------------------
if isstruct(ExpSpec) || ~isnumeric(ExpSpec)
  error('Second parameter must be experimental data.');
end

% Scale experimental spectrum
ExpSpecScaled = ExpSpec/max(abs(ExpSpec));


% Vary structure
%--------------------------------------------------------------------
% Make sure user provides one Vary structure for each Sys
if ~iscell(Vary), Vary = {Vary}; end
if numel(Vary)~=nSystems
  error(sprintf('%d spin systems given, but %d vary structure.\n Give %d vary structures.',nSystems,numel(Vary),nSystems));
end
for s=1:nSystems
  if ~isstruct(Vary{s}), Vary{s} = struct; end
end

% Make sure users are fitting with the logarithm of Diff or tcorr
for s = 1:nSystems
  if (isfield(Vary{s},'tcorr') && ~isfield(Vary{s},'logtcorr')) ||...
      (~isfield(Sys0{s},'logtcorr') && isfield(Vary{s},'logtcorr'))
    error('For least-squares fitting, use logtcorr instead of tcorr both in Sys and Vary.');
  end
  if (isfield(Vary{s},'Diff') && ~isfield(Vary{s},'logDiff')) ||...
      (~isfield(Sys0{s},'logDiff') && isfield(Vary{s},'logDiff'))
    error('For least-squares fitting, use logDiff instead of Diff both in Sys and Vary.');
  end
end

% Assert consistency between System0 and Vary structures
for s = 1:nSystems
  Fields = fieldnames(Vary{s});
  for k=1:numel(Fields)
    if ~isfield(Sys0{s},Fields{k})
      error(sprintf('Field %s is given in Vary, but not in Sys0. Remove from Vary or add to Sys0.',Fields{k}));
    elseif numel(Sys0{s}.(Fields{k})) < numel(Vary{s}.(Fields{k}))
      error(['Field ' Fields{k} ' has more elements in Vary than in Sys0.']);
    end
  end
end

% count parameters and save indices into parameter vector for each
% system
nParameters = 0;
xidx = 1;
for s=1:nSystems
  [dummy,dummy,Vals] = GetParameters(Vary{s});
  nParameters = nParameters + numel(Vals);
  Sys0{s}.xidx = xidx:xidx+numel(Vals)-1;
  xidx = xidx + numel(Vals);
end

if (nParameters==0)
  error('No variable parameters to fit.');
end


% Experimental parameters
%--------------------------------------------------------------------
if isfield(Exp,'nPoints')
  if Exp.nPoints~=numel(ExpSpec);
    error('Exp.nPoints is %d, but the spectral data vector is %d long.',...
      Exp.nPoints,numel(ExpSpec));
  end
else
  Exp.nPoints = numel(ExpSpec);
end

%if ~isfield(Exp,'Range') & ~isfield(Exp,'CenterSweep')
%  error('Please specify field range, either in Exp.Range or in Exp.CenterSweep.');
%end


% Fitting options
%--------------------------------------------------------------------

if ~isfield(FitOpt,'Method')
  FitOpt.Method = '';
end

% Method for local fitting
% (used by fit_montecarlo, fit_grid, fit_genetic)
if ~isfield(FitOpt,'LocalMethod'), FitOpt.LocalMethod = 'simplex'; end

switch FitOpt.LocalMethod
  case 'simplex', FitOpt.fitlocal = @fit_simplex;
  case 'levmar', FitOpt.fitlocal = @fit_levmar;
  otherwise
    error('Unknown local method ''%s'' (in FitOpt.LocalMethod).',FitOpt.LocalMethod)
end

FitOpt.MethodID = 1;
%if Exp.Harmonic~=0
  FitOpt.TargetID = 2;
%else
%  FitOpt.TargetID = 1;
%end
keywords = strread(FitOpt.Method,'%s');
for k = 1:numel(keywords)
  switch keywords{k}
    case 'simplex',    FitOpt.MethodID = 1;
    case 'levmar',     FitOpt.MethodID = 2;
    case 'montecarlo', FitOpt.MethodID = 3;
    case 'genetic',    FitOpt.MethodID = 4;
    case 'grid',       FitOpt.MethodID = 5;
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

MethodName{1} = 'Nelder/Mead simplex';
MethodName{2} = 'Levenberg/Marquardt';
MethodName{3} = 'Monte Carlo';
MethodName{4} = 'genetic algorithm';
MethodName{5} = 'grid search';

TargetName{1} = 'function as is';
TargetName{2} = 'integral';
TargetName{3} = 'double integral';
TargetName{4} = 'derivative';
TargetName{5} = 'Fourier transform';


%------------------------------------------------------

if ~isfield(FitOpt,'N'), FitOpt.N = 20000; end

if ~isfield(FitOpt,'Plot'), FitOpt.Plot = 1; end
if ~isfield(FitOpt,'PrintLevel'), FitOpt.PrintLevel = 1; end

if ~isfield(FitOpt,'TolFun'), FitOpt.TolFun = 1e-4; end
if ~isfield(FitOpt,'TolStep'), FitOpt.TolStep = 1e-6; end
if ~isfield(FitOpt,'maxTime'), FitOpt.maxTime = inf; end
if ~isfield(FitOpt,'RandomStart'), FitOpt.RandomStart = 0; end
if ~isfield(FitOpt,'Scaling') || isempty(FitOpt.Scaling), FitOpt.Scaling = 'maxabs'; end

if ~isfield(FitOpt,'GridSize'), FitOpt.GridSize = 7; end

% Internal parameters
if ~isfield(FitOpt,'PlotStretchFactor'), FitOpt.PlotStretchFactor = 0.05; end
if ~isfield(FitOpt,'maxGridPoints'), FitOpt.maxGridPoints = 1e5; end
if ~isfield(FitOpt,'maxParameters'), FitOpt.maxParameters = 30; end
if (nParameters>FitOpt.maxParameters)
  error('Cannot fit more than %d parameters simultaneously.',...
    FitOpt.maxParameters);
end

FitOpt.MethodString = ['Method: ' MethodName{FitOpt.MethodID} '.   Target: ', TargetName{FitOpt.TargetID}, '.   Scaling: ', FitOpt.Scaling];


%--------------------------------------------------------------


bestspec = zeros(1,numel(ExpSpec));

global UserCommand
UserCommand = 0;

if (FitOpt.Plot)
  hFig = figure(363);
  set(hFig,'Name','EasySpin least-squares fitting','NumberTitle','off');
  clf(hFig)
  set(hFig,'CloseRequestFcn',...
    'global UserCommand; UserCommand = 99; drawnow; delete(gcf);');
  hButtons(1) = uicontrol(hFig,'Style','pushbutton','String','Stop',...
    'Callback','global UserCommand; UserCommand = 1;',...
    'Position',[5 5 60 20]);
  if (FitOpt.MethodID==4) || (FitOpt.MethodID==5) || (FitOpt.MethodID==3)
    % button for ending local search
    hButtons(2) = uicontrol(hFig,'Style','pushbutton','String','Quit local',...
      'Callback','global UserCommand; UserCommand = 4; drawnow;',...
      'Position',[70 5 60 20],'Visible','off');
    % button for starting local search
    hButtons(3) = uicontrol(hFig,'Style','pushbutton','String','Start local',...
      'Callback','global UserCommand; UserCommand = 3; drawnow;',...
      'Position',[70 5 60 20]);
  end
  FitOpt.hButtons = hButtons;
  drawnow
end


%===================================================================
% Fitting algorithms
%===================================================================

if FitOpt.PrintLevel
  disp('-- esfit ------------------------------------------------');
  fprintf('Simulation function:      %s\n',char(SimFunctionName));
  fprintf('Number of components:     %d\n',nSystems);
  fprintf('Number of parameters:     %d\n',nParameters);
  fprintf('Minimization method:      %s\n',MethodName{FitOpt.MethodID});
  fprintf('Residuals computed from:  %s\n',TargetName{FitOpt.TargetID});
  fprintf('Scaling mode:             %s\n',FitOpt.Scaling);
  %disp('---------------------------------------------------------');
end

if (FitOpt.RandomStart)
  startx = 2*rand(nParameters,1) - 1;
else
  startx = zeros(nParameters,1);
end

%rand('state',sum(clock*100)); % obsoleted in R2011a
switch FitOpt.MethodID
  case 1 % Nelder/Mead simplex
    bestx = fit_simplex(@rms_,startx,FitOpt,ExpSpecScaled,Sys0,Vary,Exp,SimOpt,FitOpt);
  case 2 % Levenberg/Marquardt
    FitOpt.Gradient = FitOpt.TolFun;
    bestx = fit_levmar(@residuals_,startx,FitOpt,ExpSpecScaled,Sys0,Vary,Exp,SimOpt,FitOpt);
  case 3 % Monte Carlo
    bestx = fit_montecarlo(@assess,nParameters,FitOpt,ExpSpecScaled,Sys0,Vary,Exp,SimOpt,FitOpt);
  case 4 % Evolutionary search
    bestx = fit_genetic(@assess,nParameters,FitOpt,ExpSpecScaled,Sys0,Vary,Exp,SimOpt,FitOpt);
  case 5 % Grid search
    bestx = fit_grid(@assess,nParameters,FitOpt,ExpSpecScaled,Sys0,Vary,Exp,SimOpt,FitOpt);
  otherwise
    error('Unknown optimization method in FitOpt.Method.');
end

%===================================================================
% Final stage: finish
%===================================================================

BestSpec= 0;
for iSys=1:numel(Sys0)
  BestSys{iSys} = GetSystem(Sys0{iSys},Vary{iSys},bestx(Sys0{iSys}.xidx));
end
for iSys=1:numel(Sys0)
  % Simulate spectrum, Sys.weight is taken into account by
  % the simulation function
  if isfield(BestSys{iSys},'fcn')
    [x,b_] = feval(SimFcnName,BestSys{iSys}.fcn(BestSys{iSys}),Exp,SimOpt);
  else
    [x,b_] = feval(SimFcnName,BestSys{iSys},Exp,SimOpt);
  end
  BestSpec = BestSpec + b_;
end
BestSpecScaled = rescale(BestSpec,ExpSpecScaled,FitOpt.Scaling);
BestSpec = rescale(BestSpec,ExpSpec,FitOpt.Scaling);

[rms,Residuals_] = getresiduals(BestSpecScaled(:),ExpSpecScaled(:),FitOpt.TargetID);
RmsResidual = sqrt(mean(Residuals_.^2));
MeanResidual = mean(abs(Residuals_));
MaxResidual = max(abs(Residuals_));
StdResidual = std(abs(Residuals_));
Residuals.Rms = RmsResidual;
Residuals.Mean = MeanResidual;
Residuals.Max= MaxResidual;
Residuals.Std = StdResidual;

if FitOpt.PrintLevel && (UserCommand~=99)
  disp('---------------------------------------------------------');
  disp('Best-fit parameters:');
  str = bestfitlist(BestSys,Vary);
  fprintf(str);
  fprintf('Residuals of best fit:\n    rms  %g\n    mean abs %g,  max abs  %f,  std abs  %f\n',...
    RmsResidual,MeanResidual,MaxResidual,StdResidual);
  disp('=========================================================');
end

if (FitOpt.Plot) && (UserCommand~=99)
  figure(hFig);
  if ishandle(hButtons), delete(hButtons); end
  subplot(4,1,4);
  plot(x,BestSpec(:)-ExpSpec(:));
  xlabel('magnetic field [mT]');
  h = legend('best fit - data');
  set(h,'FontSize',8);
  axis tight
  subplot(4,1,[1 2 3]);
  h=plot(x,ExpSpec,'k.-',x,BestSpec,'g');
  set(h(2),'Color',[0 0.8 0]);
  %xlabel('magnetic field [mT]');
  h = legend('data','best fit');
  set(h,'FontSize',8);  
  axis tight
  yl = ylim;
  ylim(yl+[-1 1]*diff(yl)*FitOpt.PlotStretchFactor);
end

clear global UserCommand SimFcnName smallerr bestspec

for k=1:nSystems
  if isfield(BestSys{k},'xidx')
    BestSys{k} = rmfield(BestSys{k},'xidx');
  end
end


if (nSystems==1), BestSys = BestSys{1}; end

switch (nargout)
  case 1
    varargout = {BestSys};
  case 2
    varargout = {BestSys,BestSpec};
  case 3
    varargout = {BestSys,BestSpec,Residuals};
  case 0
    varargout = {BestSys};
end

return
%===================================================================
%===================================================================
%===================================================================

function resi = residuals_(x,ExpSpec,Sys0,Vary,Exp,SimOpt,FitOpt)
[rms,resi] = assess(x,ExpSpec,Sys0,Vary,Exp,SimOpt,FitOpt);

function rms = rms_(x,ExpSpec,Sys0,Vary,Exp,SimOpt,FitOpt)
rms = assess(x,ExpSpec,Sys0,Vary,Exp,SimOpt,FitOpt);

%==========================================================================
function varargout = assess(x,ExpSpec,Sys0,Vary,Exp,SimOpt,FitOpt)

global UserCommand SimFcnName hAxes smallerr bestspec hText
persistent BestSys;

if isempty(smallerr), smallerr = inf; end

nSystems = numel(Sys0);

% Simulate spectrum ------------------------------------------
for s = 1:nSystems
  SimSystems{s} = GetSystem(Sys0{s},Vary{s},x(Sys0{s}.xidx));
end
simspec = 0;
for s = 1:nSystems
  % SimSystems{s}.weight is taken into account in the simulation function
  if isfield(SimSystems{s},'fcn')
    simspec = simspec + ...
      feval(SimFcnName,SimSystems{s}.fcn(SimSystems{s}),Exp,SimOpt);
  else
    simspec = simspec + ...
      feval(SimFcnName,SimSystems{s},Exp,SimOpt);
  end
end

% Scale simulated spectrum to experimental spectrum ----------
simspec = rescale(simspec,ExpSpec,FitOpt.Scaling);

% Compute residuals ------------------------------
[rms,residuals] = getresiduals(simspec(:),ExpSpec(:),FitOpt.TargetID);

isNewBest = rms<smallerr;

if isNewBest
  smallerr = rms;
  bestspec = simspec;
  BestSys = SimSystems;
end

if (FitOpt.Plot) && (UserCommand~=99)
  if isempty(hAxes)
    N = numel(ExpSpec);
    set(gca,'Fontsize',8);
    hAxes = plot(1:N,ExpSpec,'k.-',1:N,bestspec,'g',1:N,simspec,'r');
    %hLegend = legend('data','best fit','current');
    %set(hLegend,'FontSize',8);
    set(gca,'XTick',[]);
    maxy = max(ExpSpec);
    miny = min(ExpSpec);
    dy = maxy-miny;
    axis tight
    set(gca,'YLim',[miny maxy] + [-1 1]*FitOpt.PlotStretchFactor*dy);
    text(min(xlim),max(ylim),' ','Vert','bottom','Tag','logLine','FontSize',7);
    text(min(xlim),min(ylim),[' ' FitOpt.MethodString],'Vert','top','Tag','methodLine','FontSize',8);
    drawnow
  else
    if ishandle(hAxes)
      set(hAxes(1),'YData',ExpSpec);
      set(hAxes(2),'YData',bestspec);
      set(hAxes(3),'YData',simspec);
      %axis tight
      drawnow
    end
  end
  if isNewBest && (UserCommand~=99)
    str = bestfitlist(BestSys,Vary);
    str = [sprintf('    RMS error: %g\n',smallerr) str];
    if (UserCommand~=99)
      if ishandle(hText)
        set(hText,'String',str,'Interpreter','none');
      else
        xl = xlim; yl = ylim;
        hText = text(xl(1),yl(2),str,'Vertical','top','FontSize',8);
        set(hText,'Color',[0 0.7 0]);
      end
    end
  end
end


if (UserCommand==2)
  UserCommand = 0;
  str = bestfitlist(BestSys,Vary);
  disp('--- current best fit parameters -------------')
  fprintf(str);
  disp('---------------------------------------------')
end

out = {rms,residuals,simspec};
varargout = out(1:nargout);
return
%==========================================================================



%==========================================================================
function Sys = GetSystem(Sys0,Vary,x)

persistent Fields Indices Vals Vary0;

Compute = isempty(Fields) | ~isequal(Vary,Vary0);
if Compute
  [Fields,Indices,Vals] = GetParameters(Vary);
  Vary0 = Vary;
end

Sys = Sys0;
Shifts = x(:).*Vals;
for k = 1:numel(x)
  f = Sys.(Fields{k});
  idx = Indices(k,:);
  f(idx(1),idx(2)) = f(idx(1),idx(2)) + Shifts(k);
  Sys.(Fields{k}) = f;
end

return
%==========================================================================



%==========================================================================
function [Fields,Indices,Values] = GetParameters(Vary)
Fields = [];
Indices = [];
Values = [];
if isempty(Vary)
  return;
end
AllFields = fieldnames(Vary);
p = 1;
for iField = 1:numel(AllFields)
  FieldValue = Vary.(AllFields{iField});
  [idx1 idx2] = find(FieldValue);
  for i = 1:numel(idx1)
    Fields{p} = AllFields{iField};
    Indices(p,:) = [idx1(i) idx2(i)];
    Values(p) = FieldValue(idx1(i),idx2(i));
    p = p + 1;
  end
end
Values = Values(:);
return
%==========================================================================


%==========================================================================
% Print from Sys values of field elements that are nonzero in Vary.
function str = bestfitlist(Sys,Vary)
nSystems = numel(Sys);
str = [];
for s=1:nSystems
  AllFields = fieldnames(Vary{s});
  if numel(AllFields)==0, continue; end
  p = 1;
  for iField = 1:numel(AllFields)
    fieldname = AllFields{iField};
    FieldValue = Sys{s}.(fieldname);
    [idx1,idx2] = find(Vary{s}.(fieldname));
    idx = sortrows([idx1(:) idx2(:)]);
    singletonDims_ = sum(size(FieldValue)==1);
    for i = 1:numel(idx1)
      Fields{p} = fieldname;
      Indices(p,:) = idx(i,:);
      singletonDims(p) = singletonDims_;
      Values(p) = FieldValue(idx(i,1),idx(i,2));
      p = p + 1;
    end
  end
  if(nSystems>1)
    str = [str sprintf('    Component %d\n',s)];
  end
  nValues = p-1;
  for p = 1:nValues
    if singletonDims(p)==2
      str = [str sprintf('       %7s:   %0.7g\n',Fields{p},Values(p))];
    elseif singletonDims(p)==1
      str = [str sprintf('    %7s(%d):   %0.7g\n',Fields{p},max(Indices(p,:)),Values(p))];
    else
      str = [str sprintf('  %7s(%d,%d):   %0.7g\n',Fields{p},Indices(p,1),Indices(p,2),Values(p))];
    end
  end
  if (nargout==0)
    fprintf(str);
  end
end
return
%==========================================================================


function [rms,residuals] = getresiduals(A,B,mode)
nanA = isnan(A);
nanB = isnan(B);
nanAB = nanA | nanB;
switch mode
  case 1
    A(nanAB) = 0;
    B(nanAB) = 0;
    residuals = A - B;
  case 2
    A(nanAB) = 0;
    B(nanAB) = 0;
    residuals = cumsum(A) - cumsum(B);
  case 3
    A(nanAB) = 0;
    B(nanAB) = 0;
    residuals = cumsum(cumsum(A)) - cumsum(cumsum(B));
  case 4
    if nanAB
      error('Cannot use FFT target with dat containin NaN.');
    end
    residuals = abs(fft(A) - fft(B));
  case 5
    A(nanAB) = 0;
    B(nanAB) = 0;
    residuals = deriv(A) - deriv(B);
end
rms = sqrt(mean(residuals.^2));
