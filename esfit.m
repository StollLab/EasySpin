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
%          -target function: 'fcn', 'int', 'dint', 'diff', 'fft'
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

global smallestError errorlist
global FitData FitOpts

smallestError = [];
errorlist = [];

% Simulation function name
%--------------------------------------------------------------------
if ~ischar(SimFunctionName)
  error('First parameter must be simulation function name.');
end
FitData.SimFcnName = SimFunctionName;
FitData.lastSetID = 0;

% System structure
%--------------------------------------------------------------------
if ~iscell(Sys0), Sys0 = {Sys0}; end
nSystems = numel(Sys0);
for s = 1:nSystems
  if ~isfield(Sys0{s},'weight'), Sys0{s}.weight = 1; end
end
FitData.nSystems = nSystems;

% Experimental spectrum
%--------------------------------------------------------------------
if isstruct(ExpSpec) || ~isnumeric(ExpSpec)
  error('Second parameter must be experimental data.');
end
FitData.nSpectra = 1;
FitData.ExpSpec = ExpSpec;
FitData.ExpSpecScaled = rescale(ExpSpec,'maxabs');

% Vary structure
%--------------------------------------------------------------------
% Make sure user provides one Vary structure for each Sys
if ~iscell(Vary), Vary = {Vary}; end
if numel(Vary)~=nSystems
  error(sprintf('%d spin systems given, but %d vary structure.\n Give %d vary structures.',nSystems,numel(Vary),nSystems));
end
for iSys=1:nSystems
  if ~isstruct(Vary{iSys}), Vary{iSys} = struct; end
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

% count parameters and save indices into parameter vector for each system
FitData.nParameters = 0;
xidx = 1;
for s=1:nSystems
  [dummy,dummy,Vals] = getParameters(Vary{s});
  FitData.nParameters = FitData.nParameters + numel(Vals);
  FitData.xidx(s) = xidx;
  xidx = xidx + numel(Vals);
end
FitData.xidx(s+1) = xidx;

if (FitData.nParameters==0)
  error('No variable parameters to fit.');
end

FitData.Vary = Vary;

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

% For field sweeps, require manual field range
%if strcmp(SimFunctionName,'pepper') || strcmp(SimFunctionName,'garlic')
%  if ~isfield(Exp,'Range') && ~isfield(Exp,'CenterSweep')
%    error('Please specify field range, either in Exp.Range or in Exp.CenterSweep.');
%  end
%end

FitData.Exp = Exp;


% Fitting options
%======================================================================

if ~isfield(FitOpt,'Scaling'), FitOpt.Scaling = 'lsq0'; end

if ~isfield(FitOpt,'Method'), FitOpt.Method = ''; end
FitOpt.MethodID = 1; % simplex
FitOpt.TargetID = 1; % function as is
if isfield(Exp,'Harmonic') && (Exp.Harmonic>0)
  FitOpt.TargetID = 2; % integral
else
  if strcmp(SimFunctionName,'pepper') || strcmp(SimFunctionName,'garlic')
    FitOpt.TargetID = 2; % integral
  end
end

keywords = strread(FitOpt.Method,'%s');
for k = 1:numel(keywords)
  switch keywords{k}
    case 'simplex',    FitOpt.MethodID = 1;
    case 'levmar',     FitOpt.MethodID = 2;
    case 'montecarlo', FitOpt.MethodID = 3;
    case 'genetic',    FitOpt.MethodID = 4;
    case 'grid',       FitOpt.MethodID = 5;
    case 'swarm',      FitOpt.MethodID = 6;
      
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

MethodNames{1} = 'Nelder/Mead simplex';
MethodNames{2} = 'Levenberg/Marquardt';
MethodNames{3} = 'Monte Carlo';
MethodNames{4} = 'genetic algorithm';
MethodNames{5} = 'grid search';
MethodNames{6} = 'particle swarm';
FitData.MethodNames = MethodNames;

TargetNames{1} = 'data as is';
TargetNames{2} = 'integral';
TargetNames{3} = 'double integral';
TargetNames{4} = 'derivative';
TargetNames{5} = 'Fourier transform';
FitData.TargetNames = TargetNames;

ScalingNames{1} = 'scale & shift (min/max)';
ScalingNames{2} = 'scale only (max abs)';
ScalingNames{3} = 'scale only (lsq)';
ScalingNames{4} = 'scale & shift (lsq0)';
ScalingNames{5} = 'scale & linear baseline (lsq1)';
ScalingNames{6} = 'scale & quad. baseline (lsq2)';
FitData.ScalingNames = ScalingNames;

ScalingString{1} = 'minmax';
ScalingString{2} = 'maxabs';
ScalingString{3} = 'lsq';
ScalingString{4} = 'lsq0';
ScalingString{5} = 'lsq1';
ScalingString{6} = 'lsq2';
FitData.ScalingString = ScalingString;

FitOpt.ScalingID = find(strcmp(FitOpt.Scaling,ScalingString));
if isempty(FitOpt.ScalingID)
  error('Unknown ''%s'' in FitOpt.Scaling.',FitOpt.Scaling);
end

%------------------------------------------------------
if ~isfield(FitOpt,'Plot'), FitOpt.Plot = 1; end
if (nargout>0), FitData.GUI = 0; else FitData.GUI = 1; end

if ~isfield(FitOpt,'PrintLevel'), FitOpt.PrintLevel = 1; end

if ~isfield(FitOpt,'nTrials'), FitOpt.nTrials = 20000; end

if ~isfield(FitOpt,'TolFun'), FitOpt.TolFun = 1e-4; end
if ~isfield(FitOpt,'TolStep'), FitOpt.TolStep = 1e-6; end
if ~isfield(FitOpt,'maxTime'), FitOpt.maxTime = inf; end
if ~isfield(FitOpt,'RandomStart'), FitOpt.RandomStart = 0; end

if ~isfield(FitOpt,'GridSize'), FitOpt.GridSize = 7; end

% Internal parameters
if ~isfield(FitOpt,'PlotStretchFactor'), FitOpt.PlotStretchFactor = 0.05; end
if ~isfield(FitOpt,'maxGridPoints'), FitOpt.maxGridPoints = 1e5; end
if ~isfield(FitOpt,'maxParameters'), FitOpt.maxParameters = 30; end
if (FitData.nParameters>FitOpt.maxParameters)
  error('Cannot fit more than %d parameters simultaneously.',...
    FitOpt.maxParameters);
end
FitData.activeParams = logical(ones(1,FitData.nParameters));

FitData.Sys0 = Sys0;
FitData.SimOpt = SimOpt;
FitOpt.IterationPrintFunction = @iterationprint;
FitOpts = FitOpt;

%=====================================================================
% Setup UI
%=====================================================================
if (FitData.GUI)
  
  % main figure
  %------------------------------------------------------------------
  hFig = findobj('Tag','esfitFigure');
  if isempty(hFig)
    hFig = figure('Tag','esfitFigure','WindowStyle','normal');
  else
    clf(hFig);
  end
  
  sz = [1000 600]; % figure size
  screensize = get(0,'ScreenSize');
  xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the screen horizontally
  ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the screen vertically
  set(hFig,'position',[xpos, ypos, sz(1), sz(2)],'units','pixels');
  set(hFig,'WindowStyle','normal','DockControls','off','MenuBar','none');
  set(hFig,'Resize','off');
  set(hFig,'Name','EasySpin least-squares fitting','NumberTitle','off');
  set(hFig,'CloseRequestFcn',...
    'global UserCommand; UserCommand = 99; drawnow; delete(gcf);');
  
  % axes
  %-----------------------------------------------------------------
  % data display
  hAx = axes('Units','pixels','Position',[50 30 600 560],'FontSize',8);
  h = plot(hAx,1,NaN,'k.-',1,NaN,'g',1,NaN,'r');
  dispData = FitData.ExpSpecScaled;
  set(h(1),'Tag','expdata','XData',1:numel(dispData),'YData',dispData);
  set(h(2),'Tag','bestsimdata','Color',[0 0.8 0]);
  set(h(3),'Tag','currsimdata');
  maxy = max(dispData); miny = min(dispData);
  set(hAx,'XLim',[1 numel(dispData)]);
  set(hAx,'YLim',[miny maxy] + [-1 1]*FitOpt.PlotStretchFactor*(maxy-miny));
  set(hAx,'XTick',[],'YTick',[]);
  
  text(mean(xlim),mean(ylim),'Press "Start" to start.',...
    'Parent',hAx,'Tag','Welcome','Margin',20,...
    'BackgroundColor',[0.7 0.9 0.7],'VerticalA','middle','HorizontalA','center');
  
  % error display
  axes('Units','pixels','Position',[670 490 100 80],'FontSize',6,'Layer','top');
  h = plot(1,NaN,'.');
  set(h,'Tag','errorline',...
    'MarkerSize',5,'Color',[0.2 0.2 0.8]);
  set(gca,'YScale','lin','XTickLabel','','YAxisLoc','right','Layer','top');
  title('log10(error)','Color','k');
  
  % axes controls
  %-----------------------------------------------------------------
  uicontrol('Style','pushbutton',...
    'Tag','lineStyleButton',...
    'Position',[50 7 20 20],...
    'String','.-',...
    'Callback',@lineStyleButtonCallback);
  
  % iteration and rms error displays
  %-----------------------------------------------------------------
  hIterText = uicontrol('Style','text','Position',[670 450 250 20]);
  set(hIterText,'FontSize',7,'Tag','logLine');
  set(hIterText,'Horizontal','left','BackgroundColor',get(gcf,'Color'));
  
  hRmsText = uicontrol('Style','text','Position',[660 330 130 30]);
  set(hRmsText,'FontSize',8,'String','','ForegroundColor',[0 0.7 0]);
  set(hRmsText,'Tag','RmsText','HorizontalAl','left');
  
  % parameter table
  %-----------------------------------------------------------------
  columnname = {'Vary','Name','best','curr'};
  columnformat = {'logical','char','char','char'};
  colEditable = [true false false false];
  [Fields,Indices,Vals] = getParameters(Vary{1});
  for f = 1:numel(Fields)
    data{f,1} = true;
    data{f,2} = [Fields{f} '(' sprintf('%d,%d',Indices(f,1),Indices(f,2)) ')'];
    %data{f,3} = sprintf('%0.6g',Vals(f));
    data{f,4} = '-';
    data{f,3} = sprintf('%0.6g',Vals(f));
  end
  %@@
  uitable('Tag','ParameterTable',...
    'Position',[660 150 330 150],...
    'ColumnFormat',columnformat,...
    'ColumnName',columnname,...
    'ColumnEditable',colEditable,...
    'ColumnWidth',{30,'auto','auto','auto'},...
    'RowName',[],...
    'Data',data);

  % set list
  %-----------------------------------------------------------------
  uicontrol(hFig,'Style','listbox',...
    'Tag','SetListBox',...
    'Position',[820 350 160 100],...
    'String','',...
    'KeyPressFcn',@deleteSetListKeyPressFcn,...
    'Tooltip','');
  uicontrol(hFig,'Style','pushbutton',...
    'Tag','deleteSetButton',...
    'Position',[820 325 40 20],...
    'String','rmv',...
    'Tooltip','remove fit set',...
    'Enable','off',...
    'Callback',@deleteSetButtonCallback);


  % popup menus
  %-----------------------------------------------------------------
  x0 = 660;
  dx = 50;
  y0 = 65; dy = 24;
  uicontrol(hFig,'Style','text',...
    'String','Method',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0+2*dy-4 dx 20]);
  uicontrol(hFig,'Style','popupmenu',...
    'Tag','MethodMenu',...
    'String',MethodNames,...
    'Value',FitOpt.MethodID,...
    'BackgroundColor','w',...
    'Tooltip','Fitting algorithm',...
    'Position',[x0+dx y0+2*dy 150 20]);
  uicontrol(hFig,'Style','text',...
    'String','Target',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0+dy-4 dx 20]);
  uicontrol(hFig,'Style','popupmenu',...
    'Tag','TargetMenu',...
    'String',TargetNames,...
    'Value',FitOpt.TargetID,...
    'BackgroundColor','w',...
    'Tooltip','Target function',...
    'Position',[x0+dx y0+dy 150 20]);
  uicontrol(hFig,'Style','text',...
    'String','Scaling',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0-4 dx 20]);
  uicontrol(hFig,'Style','popupmenu',...
    'Tag','ScalingMenu',...
    'String',ScalingNames,...
    'Value',FitOpt.ScalingID,...
    'BackgroundColor','w',...
    'Tooltip','Scaling mode',...
    'Position',[x0+dx y0 150 20]);
  
  % buttons
  %-----------------------------------------------------------------
  hei = 40;
  uicontrol(hFig,'Style','pushbutton',...
    'Tag','StartButton',...
    'String','Start',...
    'Callback',@runFitting,...
    'Visible','on',...
    'Tooltip','Start fitting',...
    'Position',[660 5 60 hei]);
  uicontrol(hFig,'Style','pushbutton',...
    'Tag','StopButton',...
    'String','Stop',...
    'Visible','off',...
    'Tooltip','Stop fitting',...
    'Callback','global UserCommand; UserCommand = 1;',...
    'Position',[660 5 60 hei]);
  drawnow
  
end

% Run fitting routine
%------------------------------------------------------------
if (~FitData.GUI)
  runFitting;
end

% Arrange outputs
%------------------------------------------------------------
if (~FitData.GUI)
  if (nSystems==1), FinalSys = FinalSys{1}; end
  switch (nargout)
    case 0, varargout = {FinalSys};
    case 1, varargout = {FinalSys};
    case 2, varargout = {FinalSys,BestSpec};
    case 3, varargout = {FinalSys,BestSpec,Residuals};
  end
else
  varargout = cell(1,nargout);
end

clear global UserCommand smallestError

%===================================================================
%===================================================================
%===================================================================

function runFitting(object,src,event)

global FitOpts FitData UserCommand
UserCommand = 0;

%===================================================================
% Update UI, pull settings from UI
%===================================================================
if (FitData.GUI)
  % Hide welcome text
  delete(findobj('Tag','Welcome'));
    
  % Hide Start button, show Stop button
  set(findobj('Tag','StopButton'),'Visible','on');
  set(findobj('Tag','StartButton'),'Visible','off');
  
  % Disable listboxes and parameter table
  set(findobj('Tag','MethodMenu'),'Enable','off');
  set(findobj('Tag','TargetMenu'),'Enable','off');
  set(findobj('Tag','ScalingMenu'),'Enable','off');
  set(findobj('Tag','ParameterTable'),'ColumnEditable',[false false false false]);
  
  % Determine selected method, target, and scaling
  FitOpts.MethodID = get(findobj('Tag','MethodMenu'),'Value');
  FitOpts.TargetID = get(findobj('Tag','TargetMenu'),'Value');
  FitOpts.Scaling = FitData.ScalingString{get(findobj('Tag','ScalingMenu'),'Value')};
  
end

%===================================================================
% Run fitting algorithm
%===================================================================

if (~FitData.GUI)
  if FitOpts.PrintLevel
    disp('-- esfit ------------------------------------------------');
    fprintf('Simulation function:      %s\n',FitData.SimFcnName);
    fprintf('Number of spectra:        %d\n',FitData.nSpectra);
    fprintf('Number of components:     %d\n',FitData.nSystems);
    fprintf('Number of parameters:     %d\n',FitData.nParameters);
    fprintf('Minimization method:      %s\n',FitData.MethodNames{FitOpts.MethodID});
    fprintf('Residuals computed from:  %s\n',FitData.TargetNames{FitOpts.TargetID});
    fprintf('Scaling mode:             %s\n',FitOpts.Scaling);
    disp('---------------------------------------------------------');
  end
end

FitData.bestspec = ones(1,numel(FitData.ExpSpec))*NaN;

if (FitOpts.RandomStart)
  startx = 2*rand(FitData.nParameters,1) - 1;
else
  startx = zeros(FitData.nParameters,1);
end

if (FitData.GUI)
  data = get(findobj('Tag','ParameterTable'),'Data');
  for iPar = 1:FitData.nParameters
    FitData.activeParams(iPar) = data{iPar,1}==1;
  end
end

startx(~FitData.activeParams) = [];
nParameters = numel(startx);

if (nParameters>0)
  switch FitOpts.MethodID
    case 1 % Nelder/Mead simplex
      bestx = fit_simplex(@assess,startx,FitOpts,FitData.ExpSpecScaled,FitData.Sys0,FitData.Vary,FitData.Exp,FitData.SimOpt,FitOpts);
    case 2 % Levenberg/Marquardt
      FitOpts.Gradient = FitOpts.TolFun;
      bestx = fit_levmar(@residuals_,startx,FitOpts,FitData.ExpSpecScaled,FitData.Sys0,FitData.Vary,FitData.Exp,FitData.SimOpt,FitOpts);
    case 3 % Monte Carlo
      bestx = fit_montecarlo(@assess,nParameters,FitOpts,FitData.ExpSpecScaled,FitData.Sys0,FitData.Vary,FitData.Exp,FitData.SimOpt,FitOpts);
    case 4 % Genetic
      bestx = fit_genetic(@assess,nParameters,FitOpts,FitData.ExpSpecScaled,FitData.Sys0,FitData.Vary,FitData.Exp,FitData.SimOpt,FitOpts);
    case 5 % Grid search
      bestx = fit_grid(@assess,nParameters,FitOpts,FitData.ExpSpecScaled,FitData.Sys0,FitData.Vary,FitData.Exp,FitData.SimOpt,FitOpts);
    case 6 % Particle swarm
      bestx = fit_particleswarm(@assess,nParameters,FitOpts,FitData.ExpSpecScaled,FitData.Sys0,FitData.Vary,FitData.Exp,FitData.SimOpt,FitOpts);
  end
else
  bestx = startx;
end

%===================================================================
% Final stage: finish
%===================================================================

BestSpec = 0;
for iSys=1:numel(FitData.Sys0)
  paramidx = FitData.xidx(iSys):FitData.xidx(iSys+1)-1;
  inactive = FitData.inactiveParams(paramidx);
  xx = zeros(1,numel(paramidx));
  xx(~inactive) = bestx;
  FinalSys{iSys} = getSystem(FitData.Sys0{iSys},FitData.Vary{iSys},xx(paramidx));
end
for iSys=1:numel(FitData.Sys0)
  % Simulate spectrum, Sys.weight is taken into account by the simulation function
  if isfield(FinalSys{iSys},'fcn')
    [x,b_] = feval(FitData.SimFcnName,FinalSys{iSys}.fcn(FinalSys{iSys}),Exp,SimOpt);
  else
    [x,b_] = feval(FitData.SimFcnName,FinalSys{iSys},FitData.Exp,FitData.SimOpt);
  end
  BestSpec = BestSpec + b_;
end
BestSpecScaled = rescale(BestSpec,FitData.ExpSpecScaled,FitOpts.Scaling);
BestSpec = rescale(BestSpec,FitData.ExpSpec,FitOpts.Scaling);

Residuals_ = getResiduals(BestSpecScaled(:),FitData.ExpSpecScaled(:),FitOpts.TargetID);
RmsResidual = sqrt(mean(Residuals_.^2));
Residuals.Rms = RmsResidual;


% Output
%===================================================================
if (~FitData.GUI)
  if FitOpts.PrintLevel && (UserCommand~=99)
    disp('---------------------------------------------------------');
    disp('Best-fit parameters:');
    str = bestfitlist(FinalSys,FitData.Vary);
    fprintf(str);
    fprintf('Residuals of best fit:\n    rms  %g\n',RmsResidual);
    disp('=========================================================');
  end
end

if FitData.GUI
  
  % Add current set to set list
  ID = FitData.lastSetID+1;
  FitData.lastSetID = ID;
  newFitSet.ID = ID;
  newFitSet.rmsd = RmsResidual;
  newFitSet.BestSpec = BestSpecScaled;
  newFitSet.Residuals = Residuals_;
  
  FitData.FitSets(ID) = newFitSet;
  
  h = findobj('Tag','SetListBox');
  str = get(h,'String');
  str{end+1} = sprintf('%d. rmsd %f',newFitSet.ID,newFitSet.rmsd);
  set(h,'String',str);
  set(h,'Value',length(str));
  set(findobj('Tag','deleteSetButton'),'Enable','on');

  % Remove current values from parameter table
  h = findobj('Tag','ParameterTable');
  Data = get(h,'Data');
  for p = 1:size(Data,1), Data{p,4} = '-'; end
  set(h,'Data',Data);
  
  % Hide current sim plot in data axes
  set(findobj('Tag','currsimdata'),'YData',NaN*ones(1,numel(BestSpec)));

  % Reactivate UI components
  set(findobj('Tag','StopButton'),'Visible','off');
  set(findobj('Tag','StartButton'),'Visible','on');
  set(findobj('Tag','MethodMenu'),'Enable','on');
  set(findobj('Tag','TargetMenu'),'Enable','on');
  set(findobj('Tag','ScalingMenu'),'Enable','on');
  set(findobj('Tag','ParameterTable'),'ColumnEditable',[true false false false]);
end

return
%===================================================================
%===================================================================
%===================================================================

function resi = residuals_(x,ExpSpec,Sys0,Vary,Exp,SimOpt,FitOpt)
[rms,resi] = assess(x,ExpSpec,Sys0,Vary,Exp,SimOpt,FitOpt);

%==========================================================================
function varargout = assess(x,ExpSpec,Sys0,Vary,Exp,SimOpt,FitOpt)

global UserCommand FitData smallestError errorlist
persistent BestSys;

if isempty(smallestError), smallestError = inf; end

% Simulate spectra ------------------------------------------
simspec = 0;
for s = 1:numel(Sys0)
  paramidx = FitData.xidx(s):FitData.xidx(s+1)-1;
  inactive = FitData.inactiveParams(paramidx);
  xx = zeros(1,numel(paramidx));
  xx(~inactive) = x;
  SimSystems{s} = getSystem(Sys0{s},Vary{s},xx(paramidx),inactive);
  % (SimSystems{s}.weight is taken into account in the simulation function)
  if isfield(SimSystems{s},'fcn')
    newsim_ = feval(FitData.SimFcnName,SimSystems{s}.fcn(SimSystems{s}),Exp,SimOpt);
  else
    newsim_ = feval(FitData.SimFcnName,SimSystems{s},Exp,SimOpt);
  end
  simspec = simspec + newsim_;
end

% Scale simulated spectrum to experimental spectrum ----------
simspec = rescale(simspec,ExpSpec,FitOpt.Scaling);

% Compute residuals ------------------------------
residuals = getResiduals(simspec(:),ExpSpec(:),FitOpt.TargetID);
rms = sqrt(mean(residuals.^2));

errorlist = [errorlist rms];
isNewBest = rms<smallestError;

if isNewBest
  smallestError = rms;
  FitData.bestspec = simspec;
  BestSys = SimSystems;
end

% update GUI
%-----------------------------------------------------------
if (FitData.GUI) && (UserCommand~=99)
  
  % update graph
  set(findobj('Tag','expdata'),'XData',1:numel(ExpSpec),'YData',ExpSpec);
  set(findobj('Tag','bestsimdata'),'XData',1:numel(ExpSpec),'YData',FitData.bestspec);
  set(findobj('Tag','currsimdata'),'XData',1:numel(ExpSpec),'YData',simspec);
  drawnow
  
  if (UserCommand~=99)
    
    % current system is new best
    if isNewBest
      [str,values] = bestfitlist(BestSys,Vary);
      str = [sprintf('RMS error: %g\n',smallestError)];
      hRmsText = findobj('Tag','RmsText');
      if ishandle(hRmsText)
        set(hRmsText,'String',str);
      end
      hParamTable = findobj('Tag','ParameterTable');
      if ishandle(hParamTable)
        data = get(hParamTable,'data');
        for p=1:numel(values)
          olddata = striphtml(data{p,3});
          newdata = sprintf('%0.6g',values(p));
          idx = 1;
          while (idx<=length(olddata)) && (idx<=length(newdata))
            if olddata(idx)~=newdata(idx), break; end
            idx = idx + 1;
          end
          data{p,3} = ['<html><font color="#000000">' newdata(1:idx-1) '</font><font color="#009900">' newdata(idx:end) '</font></html>'];
        end
        set(hParamTable,'Data',data);
      end
    end
    
    % current system
    [str,values] = bestfitlist(SimSystems,Vary);
    hParamTable = findobj('Tag','ParameterTable');
    if ishandle(hParamTable)
      data = get(hParamTable,'data');
      for p=1:numel(values)
        olddata = striphtml(data{p,4});
        newdata = sprintf('%0.6f',values(p));
        idx = 1;
        while (idx<=length(olddata)) && (idx<=length(newdata))
          if olddata(idx)~=newdata(idx), break; end
          idx = idx + 1;
        end
        data{p,4} = ['<html><font color="#000000">' newdata(1:idx-1) '</font><font color="#ff0000">' newdata(idx:end) '</font></html>'];
      end
      set(hParamTable,'Data',data);
    end
  end
  
  hErrorLine = findobj('Tag','errorline');
  if ~isempty(hErrorLine)
    n = min(200,numel(errorlist));
    set(hErrorLine,'XData',1:n,'YData',log10(errorlist(end-n+1:end)));
    axis(get(hErrorLine,'Parent'),'tight');
    drawnow
  end
  
end
%-------------------------------------------------------------------

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
function Sys = getSystem(Sys0,Vary,x,skip)

persistent Vary0 Fields Indices Vals;

if isempty(Vary0) || ~isequal(Vary,Vary0);
  Vary0 = Vary;
  [Fields,Indices,Vals] = getParameters(Vary);
end

Sys = Sys0;
Shifts = x(:).*Vals;
if nargin>3, Shifts(skip) = 0; end
for k = 1:numel(x)
  f = Sys.(Fields{k});
  idx = Indices(k,:);
  f(idx(1),idx(2)) = f(idx(1),idx(2)) + Shifts(k);
  Sys.(Fields{k}) = f;
end

return
%==========================================================================



%==========================================================================
function [Fields,Indices,Values] = getParameters(Vary)
Fields = [];
Indices = [];
Values = [];
if isempty(Vary), return; end
allFields = fieldnames(Vary);
p = 1;
for iField = 1:numel(allFields)
  FieldValue = Vary.(allFields{iField});
  [idx1,idx2] = find(FieldValue);
  for i = 1:numel(idx1)
    Fields{p} = allFields{iField};
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
function [str,Values] = bestfitlist(Sys,Vary)
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
  if (nSystems>1)
    str = [str sprintf('  Component %d\n',s)];
  end
  nValues = p-1;
  for p = 1:nValues
    if singletonDims(p)==2
      str = [str sprintf('     %7s:   %0.7g\n',Fields{p},Values(p))];
    elseif singletonDims(p)==1
      str = [str sprintf('  %7s(%d):   %0.7g\n',Fields{p},max(Indices(p,:)),Values(p))];
    else
      str = [str sprintf('%7s(%d,%d):   %0.7g\n',Fields{p},Indices(p,1),Indices(p,2),Values(p))];
    end
  end
  if (nargout==0)
    fprintf(str);
  end
end
return
%==========================================================================


%==========================================================================
function residuals = getResiduals(A,B,mode)
residuals = A - B;
idxNaN = isnan(A) | isnan(B);
residuals(idxNaN) = 0; % ignore NaNs in either A or B
switch mode
  case 1 % fcn
  case 2 % int
    residuals = cumsum(residuals);
  case 3 % iint
    residuals = cumsum(cumsum(residuals));
  case 4 % fft
    residuals = abs(fft(residuals));
  case 5 % diff
    residuals = deriv(residuals);
end
return
%==========================================================================

function iterationprint(str)
hLogLine = findobj('Tag','logLine');
if isempty(hLogLine)
  disp(str);
else
  set(hLogLine,'String',str);
end

function str = striphtml(str)
html = 0;
for k=1:numel(str)
  if ~html
    rmv(k) = false;
    if str(k)=='<', html = 1; rmv(k) = true; end
  else
    rmv(k) = true;
    if str(k)=='>', html = 0; end
  end
end
str(rmv) = [];
return

function plotFittingResult
if (FitOpt.Plot) && (UserCommand~=99)
  close(hFig); clf
  
  subplot(4,1,4);
  plot(x,BestSpec(:)-ExpSpec(:));
  h = legend('best fit - data');
  legend boxoff
  set(h,'FontSize',8);
  axis tight
  height4 = get(gca,'Position'); height4 = height4(4);
  
  subplot(4,1,[1 2 3]);
  h = plot(x,ExpSpec,'k.-',x,BestSpec,'g');
  set(h(2),'Color',[0 0.8 0]);
  h = legend('data','best fit');
  legend boxoff
  set(h,'FontSize',8);
  axis tight
  yl = ylim;
  yl = yl+[-1 1]*diff(yl)*FitOpt.PlotStretchFactor;
  ylim(yl);
  height123 = get(gca,'Position'); height123 = height123(4);
  
  subplot(4,1,4);
  yl = ylim;
  ylim(mean(yl)+[-1 1]*diff(yl)*height123/height4/2);
  
end
return

function lineStyleButtonCallback(object,src,event)
h = findobj('Tag','ExpData');
if get(h,'LineStyle')=='-'
  set(h,'LineStyle','none');
  set(h,'Marker','.');
else
  set(h,'LineStyle','-');
  set(h,'Marker','.');
end
return

function deleteSetButtonCallback(object,src,event)
h = findobj('Tag','SetListBox');
idx = get(h,'Value');
str = get(h,'String');
if ~isempty(str)
  str(idx) = [];
  if idx>length(str), idx = length(str); end
  if (idx==0), idx = 1; end
  set(h,'Value',idx);
  set(h,'String',str);
end
if isempty(str)
  set(findobj('Tag','deleteSetButton'),'Enable','off');
end
return

function deleteSetListKeyPressFcn(object,event)
if strcmp(event.Key,'delete')
  deleteSetButtonCallback(object,gco,event);
end
return
