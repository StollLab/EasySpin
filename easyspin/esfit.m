% esfit   least-squares fitting for EPR spectral simulations
%
%   esfit(simfunc,expspc,Sys0,Vary,Exp)
%   esfit(simfunc,expspc,Sys0,Vary,Exp,SimOpt)
%   esfit(simfunc,expspc,Sys0,Vary,Exp,SimOpt,FitOpt)
%   bestsys = esfit(...)
%   [bestsys,bestspc] = esfit(...)
%   [bestsys,bestspc,residuals] = esfit(...)
%
% Input:
%     simfunc     simulation function handle (@pepper, @garlic, @salt, ...
%                   @chili, or user-defined function)
%     expspc      experimental spectrum, a vector of data points
%     Sys0        starting values for spin system parameters
%     Vary        allowed variation of parameters
%     Exp         experimental parameter, for simulation function
%     SimOpt      options for the simulation algorithms
%     FitOpt      options for the fitting algorithms
%        Method   string containing kewords for
%          -algorithm: 'simplex','levmar','montecarlo','genetic','grid',
%                      'swarm'
%          -target function: 'fcn', 'int', 'dint', 'diff', 'fft'
%        Scaling  string with scaling method keyword
%          'maxabs' (default), 'minmax', 'lsq', 'lsq0','lsq1','lsq2','none'
%        OutArg   two numbers [nOut iOut], where nOut is the number of
%                 outputs of the simulation function and iOut is the index
%                 of the output argument to use for fitting
% Output:
%     bestsys     spin system with fitted parameters
%     bestspc     fitted simulated spectrum
%     residuals   residuals between fitted and experimental spectrum (vector)

function varargout = esfit(SimFunction,ExpSpec,Sys0,Vary,Exp,SimOpt,FitOpt)

if nargin==1 && isnumeric(SimFunction) && SimFunction==1
  SimFunction = @pepper;
  Sys0.g = [1.9 1.97 2.1];
  Sys0.lwpp = 3;
  
  Exp.mwFreq = 9.5;
  Exp.Range = [300 380];
  
  [~,spc] = pepper(Sys0,Exp);
  ExpSpec = 15*addnoise(spc,50,'n');
  
  Sys0.g = Sys0.g + (rand(1,3)-0.5)*0.02;
  
  Vary.g = [1 1 1]*0.04;
  Vary.lwpp = 0.9*Sys0.lwpp;
  
  SimOpt = struct;
  
end

if nargin==0, help(mfilename); return; end

if nargin<5
  error('Not enough input arguments. At least 5 are required.')
end
if nargin>7
  error('Too many input arguments. No more than 7 are accepted.')
end

% Check expiry date
error(eschecker);

%if nargin<5, error('Not enough inputs.'); end
if nargin<6, SimOpt = struct; end
if nargin<7, FitOpt = struct; end

if isempty(FitOpt), FitOpt = struct; end
if ~isstruct(FitOpt)
  error('FitOpt (7th input argument of esfit) must be a structure.');
end

global FitData FitOpts
FitData = [];
FitOpts = [];
FitData.currFitSet = [];

% Simulation function
%--------------------------------------------------------------------
if ~isa(SimFunction,'function_handle')
  str = 'The simulation function (1st input) must be a function handle (with @).';
  if ischar(SimFunction)
    error('%s\nUse esfit(@%s,...) instead of esfit(''%s'',...).',str,SimFunction,SimFunction);
  else
    error('%s\nFor example, to use the function pepper(...), use esfit(@pepper,...).',str);
  end
end

try
  nargin(SimFunction);
catch
  error('The function given as 1st input cannot be found.');
end

FitData.SimFcnName = func2str(SimFunction);
FitData.SimFcn = SimFunction;

FitData.lastSetID = 0;

% System structure
%-------------------------------------------------------------------------------
if ~iscell(Sys0), Sys0 = {Sys0}; end
nSystems = numel(Sys0);
for s = 1:nSystems
  if ~isfield(Sys0{s},'weight'), Sys0{s}.weight = 1; end
end
FitData.nSystems = nSystems;

% Experimental spectrum
%-------------------------------------------------------------------------------
if isstruct(ExpSpec) || ~isnumeric(ExpSpec)
  error('Second parameter must be experimental data.');
end
FitData.nSpectra = 1;
FitData.ExpSpec = ExpSpec;
FitData.ExpSpecScaled = rescaledata(ExpSpec,'maxabs');

% Vary structure
%-------------------------------------------------------------------------------
% Make sure user provides one Vary structure for each Sys
if ~iscell(Vary), Vary = {Vary}; end
if numel(Vary)~=nSystems
  error('%d spin systems given, but %d vary structure.\n Give %d vary structures.',nSystems,numel(Vary),nSystems);
end
for iSys = 1:nSystems
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
  for k = 1:numel(Fields)
    if ~isfield(Sys0{s},Fields{k})
      error('Field %s is given in Vary, but not in Sys0. Remove from Vary or add to Sys0.',Fields{k});
    elseif numel(Sys0{s}.(Fields{k})) < numel(Vary{s}.(Fields{k}))
      error(['Field ' Fields{k} ' has more elements in Vary than in Sys0.']);
    end
  end
  clear Fields
end

% Count parameters and save indices into parameter vector for each system
for iSys = 1:nSystems
  [~,~,v_] = getParameters(Vary{iSys});
  VaryVals(iSys) = numel(v_);
end
FitData.xidx = cumsum([1 VaryVals]);
FitData.nParameters = sum(VaryVals);

if FitData.nParameters==0
  error('No variable parameters to fit.');
end

FitData.Vary = Vary;

% Experimental parameters
%-------------------------------------------------------------------------------
if isfield(Exp,'nPoints')
  if Exp.nPoints~=numel(ExpSpec)
    error('Exp.nPoints is %d, but the spectral data vector is %d long.',...
      Exp.nPoints,numel(ExpSpec));
  end
else
  Exp.nPoints = numel(ExpSpec);
end

% For field and frequency sweeps, require manual field range (to prevent
% users from comparing sim and exp spectra with different ranges)
if strcmp(SimFunction,{'pepper','garlic','chili','salt'})
  if ~any(isfield(Exp,{'Range','CenterSweep','mwRange','mwCenterSweep'}))
    error('Please specify field or frequency range, in Exp.Range/Exp.mwRange or in Exp.CenterSweep/Exp.mwCenterSweep.');
  end
end

FitData.Exp = Exp;


% Fitting options
%===============================================================================
if ~isfield(FitOpt,'OutArg')
  FitData.nOutArguments = abs(nargout(FitData.SimFcn));
  FitData.OutArgument = FitData.nOutArguments;
else
  if numel(FitOpt.OutArg)~=2
    error('FitOpt.OutArg must contain two values [nOut iOut]');
  end
  if FitOpt.OutArg(2)>FitOpt.OutArg(1)
    error('FitOpt.OutArg: second number cannot be larger than first one.');
  end
  FitData.nOutArguments = FitOpt.OutArg(1);
  FitData.OutArgument = FitOpt.OutArg(2);
  
end

if ~isfield(FitOpt,'Scaling'), FitOpt.Scaling = 'lsq0'; end

if ~isfield(FitOpt,'Method'), FitOpt.Method = ''; end
FitOpt.MethodID = 1; % simplex
FitOpt.TargetID = 1; % function as is
if isfield(Exp,'Harmonic') && (Exp.Harmonic>0)
  FitOpt.TargetID = 2; % integral
else
  if strcmp(FitData.SimFcnName,'pepper') || strcmp(FitData.SimFcnName,'garlic')
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
ScalingNames{7} = 'no scaling';
FitData.ScalingNames = ScalingNames;

ScalingString{1} = 'minmax';
ScalingString{2} = 'maxabs';
ScalingString{3} = 'lsq';
ScalingString{4} = 'lsq0';
ScalingString{5} = 'lsq1';
ScalingString{6} = 'lsq2';
ScalingString{7} = 'none';
FitData.ScalingString = ScalingString;

StartpointNames{1} = 'center of range';
StartpointNames{2} = 'random within range';
StartpointNames{3} = 'selected parameter set';


FitOpt.ScalingID = find(strcmp(FitOpt.Scaling,ScalingString));
if isempty(FitOpt.ScalingID)
  error('Unknown ''%s'' in FitOpt.Scaling.',FitOpt.Scaling);
end

%-------------------------------------------------------------------------------
if ~isfield(FitOpt,'Plot'), FitOpt.Plot = true; end
if nargout>0, FitData.GUI = false; else, FitData.GUI = true; end

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
if FitData.nParameters>FitOpt.maxParameters
  error('Cannot fit more than %d parameters simultaneously.',...
    FitOpt.maxParameters);
end
FitData.inactiveParams = false(1,FitData.nParameters);

FitData.Sys0 = Sys0;
FitData.SimOpt = SimOpt;
FitOpt.IterationPrintFunction = @iterationprint;
FitOpts = FitOpt;

%===============================================================================
% Setup UI
%===============================================================================
if FitData.GUI
  clc
  
  % main figure
  %-----------------------------------------------------------------------------
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
  %-----------------------------------------------------------------
  excludedRegions = [];
  % data display
  hAx = axes('Parent',hFig,'Units','pixels',...
    'Position',[10 10 640 580],'FontSize',8,'Layer','top');
  NaNdata = ones(1,Exp.nPoints)*NaN;
  dispData = FitData.ExpSpecScaled;
  maxy = max(dispData); miny = min(dispData);
  YLimits = [miny maxy] + [-1 1]*FitOpt.PlotStretchFactor*(maxy-miny);
  for r = 1:size(excludedRegions,1)
    h = patch(excludedRegions(r,[1 2 2 1]),YLimits([1 1 2 2]),[1 1 1]*0.8);
    set(h,'EdgeColor','none');
  end
  x = 1:Exp.nPoints;
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
  %-----------------------------------------------------------------
  columnname = {'','Name','best','current','lower','upper'};
  columnformat = {'logical','char','char','char','char','char'};
  colEditable = [true false false false true true];
  [FitData.parNames,FitData.CenterVals,FitData.VaryVals] = getParamList(Sys0,Vary);
  data = cell(numel(FitData.parNames),6);
  for p = 1:numel(FitData.parNames)
    data{p,1} = true;
    data{p,2} = FitData.parNames{p};
    data{p,3} = '-';
    data{p,4} = '-';
    p0 = FitData.CenterVals(p);
    dp = FitData.VaryVals(p);
    lower = p0-dp;
    upper = p0+dp;
    data{p,5} = sprintf('%0.6g',lower);
    data{p,6} = sprintf('%0.6g',upper);
  end
  x0 = 660; y0 = 400; dx = 80;
  uitable('Tag','ParameterTable',...
    'FontSize',8,...
    'Position',[x0 y0 330 150],...
    'ColumnFormat',columnformat,...
    'ColumnName',columnname,...
    'ColumnEditable',colEditable,...
    'CellEditCallback',@tableEditCallback,...
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
  h = uicontrol(hFig,'Style','text',...
    'String','tbd',...
    'ForeGroundColor','b',...
    'Tooltip','Simulation function name',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0+dx y0+170 dx 20]);
  set(h,'Tooltip',sprintf('using output no. %d of %d',FitData.nOutArguments,FitData.OutArgument));
  set(h,'String',FitData.SimFcnName);

  % popup menus
  %-----------------------------------------------------------------------------
  x0 = 660; dx = 60; y0 = 290; dy = 24;
  uicontrol(hFig,'Style','text',...
    'String','Method',...
    'FontWeight','bold',...
    'HorizontalAlign','left',...
    'BackgroundColor',get(gcf,'Color'),...
    'Position',[x0 y0+3*dy-4 dx 20]);
  uicontrol(hFig,'Style','popupmenu',...
    'Tag','MethodMenu',...
    'String',MethodNames,...
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
    'String',TargetNames,...
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
    'String',ScalingNames,...
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
    'String',StartpointNames,...
    'Value',1,...
    'BackgroundColor','w',...
    'Tooltip','Starting point for fit',...
    'Position',[x0+dx y0 150 20]);
  if FitOpts.Startpoint==2, set(h,'Value',2); end
  
  % Start/Stop buttons
  %-----------------------------------------------------------------
  pos =  [x0+220 y0-3+30 110 67];
  pos1 = [x0+220 y0-3    110 30];
  uicontrol(hFig,'Style','pushbutton',...
    'Tag','StartButton',...
    'String','Start',...
    'Callback',@runFitting,...
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
  %-----------------------------------------------------------------
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

% Run fitting routine
%------------------------------------------------------------
if ~FitData.GUI
  [BestSys,BestSpec,Residuals] = runFitting;
end

% Arrange outputs
%------------------------------------------------------------
if ~FitData.GUI
  if nSystems==1, BestSys = BestSys{1}; end
  switch nargout
    case 0, varargout = {BestSys};
    case 1, varargout = {BestSys};
    case 2, varargout = {BestSys,BestSpec};
    case 3, varargout = {BestSys,BestSpec,Residuals};
  end
else
  varargout = cell(1,nargout);
end

clear global UserCommand

%===================================================================
%===================================================================
%===================================================================

function [FinalSys,BestSpec,Residuals] = runFitting(~,~,~)

global FitOpts FitData UserCommand

UserCommand = 0;

%===================================================================
% Update UI, pull settings from UI
%===================================================================
if FitData.GUI
    
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
  FitOpts.MethodID = get(findobj('Tag','MethodMenu'),'Value');
  FitOpts.TargetID = get(findobj('Tag','TargetMenu'),'Value');
  FitOpts.Scaling = FitData.ScalingString{get(findobj('Tag','ScalingMenu'),'Value')};
  FitOpts.Startpoint = get(findobj('Tag','StartpointMenu'),'Value');
  
end

%===================================================================
% Run fitting algorithm
%===================================================================

if ~FitData.GUI
  if FitOpts.PrintLevel
    disp('-- esfit ------------------------------------------------');
    fprintf('Simulation function:      %s\n',FitData.SimFcnName);
    fprintf('Problem size:             %d spectra, %d components, %d parameters\n',FitData.nSpectra,FitData.nSystems,FitData.nParameters);
    fprintf('Minimization method:      %s\n',FitData.MethodNames{FitOpts.MethodID});
    fprintf('Residuals computed from:  %s\n',FitData.TargetNames{FitOpts.TargetID});
    fprintf('Scaling mode:             %s\n',FitOpts.Scaling);
    disp('---------------------------------------------------------');
  end
end

%FitData.bestspec = ones(1,numel(FitData.ExpSpec))*NaN;

if FitData.GUI
  data = get(findobj('Tag','ParameterTable'),'Data');
  for iPar = 1:FitData.nParameters
    FitData.inactiveParams(iPar) = data{iPar,1}==0;
  end
end

switch FitOpts.Startpoint
  case 1 % center of range
    startx = zeros(FitData.nParameters,1);
  case 2 % random
    startx = 2*rand(FitData.nParameters,1) - 1;
    startx(FitData.inactiveParams) = 0;
  case 3 % selected parameter set
    h = findobj('Tag','SetListBox');
    s = h.String;
    if ~isempty(s)
      s = s{h.Value};
      ID = sscanf(s,'%d');
      idx = find([FitData.FitSets.ID]==ID);
      if ~isempty(idx)
        startx = FitData.FitSets(idx).bestx;
      else
        error('Could not locate selected parameter set.');
      end
    else
      startx = zeros(FitData.nParameters,1);
    end
end
FitData.startx = startx;

x0_ = startx;
x0_(FitData.inactiveParams) = [];

bestx = startx;
if strcmp(FitOpts.Scaling, 'none')
  fitspc = FitData.ExpSpec;
else
  fitspc = FitData.ExpSpecScaled;
end

funArgs = {fitspc,FitData,FitOpts};  % input args for assess and residuals_

if FitData.GUI
  h = findobj('Tag','statusText');
  if ~strcmp(h.String,'running')
    set(h,'String','running');
    drawnow
  end
end


nParameters_ = numel(x0_);
if nParameters_>0
  switch FitOpts.MethodID
    case 1 % Nelder/Mead simplex
      bestx0_ = esfit_simplex(@rmsd_,x0_,FitOpts,funArgs{:});
    case 2 % Levenberg/Marquardt
      FitOpts.Gradient = FitOpts.TolFun;
      bestx0_ = esfit_levmar(@residuals_,x0_,FitOpts,funArgs{:});
    case 3 % Monte Carlo
      bestx0_ = esfit_montecarlo(@rmsd_,nParameters_,FitOpts,funArgs{:});
    case 4 % Genetic
      bestx0_ = esfit_genetic(@rmsd_,nParameters_,FitOpts,funArgs{:});
    case 5 % Grid search
      bestx0_ = esfit_grid(@rmsd_,nParameters_,FitOpts,funArgs{:});
    case 6 % Particle swarm
      bestx0_ = esfit_swarm(@rmsd_,nParameters_,FitOpts,funArgs{:});
  end
  bestx(~FitData.inactiveParams) = bestx0_;
end

if FitData.GUI
  set(findobj('Tag','statusText'),'String','');
end

if FitData.GUI
  
  % Remove current values from parameter table
  hTable = findobj('Tag','ParameterTable');
  Data = hTable.Data;
  for p = 1:size(Data,1), Data{p,4} = '-'; end
  set(hTable,'Data',Data);
  
  % Hide current sim plot in data axes
  set(findobj('Tag','currsimdata'),'YData',NaN*ones(1,numel(FitData.ExpSpec)));
  hErrorLine = findobj('Tag','errorline');
  set(hErrorLine,'XData',1,'YData',NaN);
  axis(hErrorLine.Parent,'tight');
  drawnow
  set(findobj('Tag','logLine'),'String','');

  % Reactivate UI components
  set(findobj('Tag','SaveButton'),'Enable','on');
  
  if isfield(FitData,'FitSets') && numel(FitData.FitSets)>0
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
  
end

%===================================================================
% Final stage: finish
%===================================================================

% compile best-fit system structures
[FinalSys,bestvalues] = getSystems(FitData.Sys0,FitData.Vary,bestx);

% Simulate best-fit spectrum
if numel(FinalSys)==1
  fs = FinalSys{1};
else
  fs = FinalSys;
end
[out{1:FitData.nOutArguments}] = FitData.SimFcn(fs,FitData.Exp,FitData.SimOpt);

% (SimSystems{s}.weight is taken into account in the simulation function)
BestSpec = out{FitData.OutArgument}; % pick last output argument
BestSpecScaled = rescaledata(BestSpec,FitData.ExpSpecScaled,FitOpts.Scaling);
BestSpec = rescaledata(BestSpec,FitData.ExpSpec,FitOpts.Scaling);

Residuals = calculateResiduals(BestSpecScaled(:),FitData.ExpSpecScaled(:),FitOpts.TargetID);
Residuals = Residuals.'; % col -> row
rmsd = sqrt(mean(Residuals.^2));

% Output
%===============================================================================
if ~FitData.GUI
  
  if FitOpts.PrintLevel && UserCommand~=99
    disp('---------------------------------------------------------');
    disp('Best-fit parameters:');
    str = bestfitlist(FinalSys,FitData.Vary);
    fprintf(str);
    fprintf('Residuals of best fit:\n    rmsd  %g\n',rmsd);
    disp('=========================================================');
  end

else
  
  % Save current set to set list
  newFitSet.rmsd = rmsd;
  if strcmp(FitOpts.Scaling, 'none')
    newFitSet.fitSpec = BestSpec;
    newFitSet.expSpec = FitData.ExpSpec;
  else
    newFitSet.fitSpec = BestSpecScaled;
    newFitSet.expSpec = FitData.ExpSpecScaled;
  end
  newFitSet.residuals = Residuals;
  newFitSet.bestx = bestx;
  newFitSet.bestvalues = bestvalues;
  TargetKey = {'fcn','int','iint','diff','fft'};
  newFitSet.Target = TargetKey{FitOpts.TargetID};
  newFitSet.Scaling = FitOpts.Scaling;
  if numel(FinalSys)==1
    newFitSet.Sys = FinalSys{1};
  else
    newFitSet.Sys = FinalSys;
  end
  FitData.currFitSet = newFitSet;
  
end

return
%===============================================================================
%===============================================================================
%===============================================================================

function resi = residuals_(x,ExpSpec,FitDat,FitOpt)
[~,resi] = rmsd_(x,ExpSpec,FitDat,FitOpt);

%===============================================================================
function varargout = rmsd_(x,ExpSpec,FitDat,FitOpt)

global UserCommand FitData
persistent BestSys;

if ~isfield(FitData,'smallestError') || isempty(FitData.smallestError)
  FitData.smallestError = inf;
end
if ~isfield(FitData,'errorlist')
  FitData.errorlist = [];
end

Sys0 = FitDat.Sys0;
Vary = FitDat.Vary;
Exp = FitDat.Exp;
SimOpt = FitDat.SimOpt;


% Simulate spectra ------------------------------------------
inactive = FitData.inactiveParams;
x_all = FitData.startx;
x_all(~inactive) = x;
[SimSystems,simvalues] = getSystems(Sys0,Vary,x_all);
try
  if numel(SimSystems)==1
    [out{1:FitData.nOutArguments}] = FitData.SimFcn(SimSystems{1},Exp,SimOpt);
  else
    [out{1:FitData.nOutArguments}] = FitData.SimFcn(SimSystems,Exp,SimOpt);
  end
catch
  % TBD: better error recovery; error display in GUI
  error(lasterr);
end
% (SimSystems{s}.weight is taken into account in the simulation function)
simspec = out{FitData.OutArgument}; % pick last output argument

% Scale simulated spectrum to experimental spectrum ----------
simspec = rescaledata(simspec,ExpSpec,FitOpt.Scaling);

% Compute residuals ------------------------------
residuals = calculateResiduals(simspec(:),ExpSpec(:),FitOpt.TargetID);
rmsd = real(sqrt(mean(residuals.^2)));

FitData.errorlist = [FitData.errorlist rmsd];
isNewBest = rmsd<FitData.smallestError;

if isNewBest
  FitData.smallestError = rmsd;
  FitData.bestspec = simspec;
  BestSys = SimSystems;
end

% update GUI
%-----------------------------------------------------------
if FitData.GUI && UserCommand~=99
  
  % update plot
  x = 1:numel(ExpSpec);
  set(findobj('Tag','expdata'),'XData',x,'YData',ExpSpec);
  set(findobj('Tag','bestsimdata'),'XData',x,'YData',real(FitData.bestspec));
  set(findobj('Tag','currsimdata'),'XData',x,'YData',real(simspec));
  
  % readjust vertical range
  dispData = [ExpSpec(:); real(FitData.bestspec(:)); real(simspec(:))];
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
    for p = 1:numel(simvalues)
      olddata = striphtml(data{p,4});
      newdata = sprintf('%0.6f',simvalues(p));
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
      [~,values] = getSystems(BestSys,Vary);
      
      str = sprintf(' best RMSD: %g\n',(FitData.smallestError));
      hRmsText = findobj('Tag','RmsText');
      set(hRmsText,'String',str);
      
      for p = 1:numel(values)
        olddata = striphtml(data{p,3});
        newdata = sprintf('%0.6g',values(p));
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
    n = min(100,numel(FitData.errorlist));
    set(hErrorLine,'XData',1:n,'YData',log10(FitData.errorlist(end-n+1:end)));
    ax = hErrorLine.Parent;
    axis(ax,'tight');
    drawnow
  end
  
end
%-------------------------------------------------------------------

if UserCommand==2
  UserCommand = 0;
  str = bestfitlist(BestSys,Vary);
  disp('--- current best fit parameters -------------')
  fprintf(str);
  disp('---------------------------------------------')
end

out = {rmsd,residuals,simspec};
varargout = out(1:nargout);
return
%==========================================================================



%==========================================================================
% Calculate spin systems with values based on Sys0 (starting points), Vary
% (parameters to vary, and their vary range), and x (current point in vary
% range)
function [Sys,values] = getSystems(Sys0,Vary,x)
global FitData
values = [];
if nargin==3, x = x(:); end
for iSys = 1:numel(Sys0)
  [Fields,Indices,VaryVals] = getParameters(Vary{iSys});
  
  if isempty(VaryVals)
    % no parameters varied in this spin system
    Sys{iSys} = Sys0{iSys};
    continue
  end
  
  thisSys = Sys0{iSys};
  
  pidx = FitData.xidx(iSys):FitData.xidx(iSys+1)-1;
  if nargin<3
    Shifts = zeros(numel(VaryVals),1);
  else
    Shifts = x(pidx).*VaryVals(:);
  end
  values_ = [];
  for p = 1:numel(VaryVals)
    f = thisSys.(Fields{p});
    idx = Indices(p,:);
    values_(p) = f(idx(1),idx(2)) + Shifts(p);
    f(idx(1),idx(2)) = values_(p);
    thisSys.(Fields{p}) = f;
  end
  
  values = [values values_];
  Sys{iSys} = thisSys;
  
end

return
%==========================================================================


%==========================================================================
function [parNames,parCenter,parVary] = getParamList(Sys,Vary)
nSystems = numel(Sys);
p = 1;
for s = 1:nSystems
  allFields = fieldnames(Vary{s});
  for iField = 1:numel(allFields)
    fieldName = allFields{iField};
    CenterValue = Sys{s}.(fieldName);
    VaryValue = Vary{s}.(fieldName);
    [idx1,idx2] = find(VaryValue);
    idx = sortrows([idx1(:) idx2(:)]);
    nValues = numel(idx1);
    for iVal = 1:nValues
      parCenter(p) = CenterValue(idx(iVal,1),idx(iVal,2));
      parVary(p) = VaryValue(idx(iVal,1),idx(iVal,2));
      Indices = idx(iVal,:);
      if isvector(CenterValue)
        idxString_ = sprintf('(%d)',max(Indices));
      elseif ismatrix(CenterValue)
        idxString_ = sprintf('(%d,%d)',Indices(1),Indices(2));
      else
        idxString_ = '';
      end
      parNames{p} = [fieldName idxString_];
      if nSystems>1
        sysID = char('A'-1+s);
        parNames{p} = [sysID '.' parNames{p}];
      end
      p = p + 1;
    end
  end
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
  Value = Vary.(allFields{iField});
  [idx1,idx2] = find(Value);
  idx = sortrows([idx1(:) idx2(:)]);
  for i = 1:numel(idx1)
    Fields{p} = allFields{iField};
    Indices(p,:) = [idx(i,1) idx(i,2)];
    Values(p) = Value(idx(i,1),idx(i,2));
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
p = 1;
for s=1:nSystems
  AllFields = fieldnames(Vary{s});
  if numel(AllFields)==0, continue; end
  for iField = 1:numel(AllFields)
    fieldname = AllFields{iField};
    FieldValue = Sys{s}.(fieldname);
    [idx1,idx2] = find(Vary{s}.(fieldname));
    idx = sortrows([idx1(:) idx2(:)]);
    singletonDims_ = sum(size(FieldValue)==1);
    for i = numel(idx1):-1:1
      Fields{p} = fieldname;
      Indices(p,:) = idx(i,:);
      singletonDims(p) = singletonDims_;
      Values(p) = FieldValue(idx(i,1),idx(i,2));
      Component(p) = s;
      p = p + 1;
    end
  end
end
nParameters = p-1;

for p = 1:nParameters
  if nSystems>1 && (p==1 || Component(p-1)~=Component(p))
    str = [str sprintf('component %s\n',char('A'-1+Component(p)))];
  end
  if singletonDims(p)==2
    str = [str sprintf('     %7s:   %0.7g\n',Fields{p},Values(p))];
  elseif singletonDims(p)==1
    str = [str sprintf('  %7s(%d):   %0.7g\n',Fields{p},max(Indices(p,:)),Values(p))];
  else
    str = [str sprintf('%7s(%d,%d):   %0.7g\n',Fields{p},Indices(p,1),Indices(p,2),Values(p))];
  end
end

if nargout==0, fprintf(str); end
return
%==========================================================================


%==========================================================================
function residuals = calculateResiduals(A,B,mode)
switch mode
  case 1 % fcn
    fcn = @(x) x;
  case 2 % int
    fcn = @(x) cumsum(x);
  case 3 % iint
    fcn = @(x) cumsum(cumsum(x));
  case 4 % fft
    fcn = @(x) abs(fft(x));
  case 5 % diff
    fcn = @(x) deriv(x);
end
residuals = fcn(A) - fcn(B);
idxNaN = isnan(A) | isnan(B);
residuals(idxNaN) = 0; % ignore NaNs in either A or B
return
%==========================================================================

%==========================================================================
function iterationprint(str)
hLogLine = findobj('Tag','logLine');
if isempty(hLogLine)
  disp(str);
else
  set(hLogLine,'String',str);
end
%==========================================================================


%==========================================================================
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
return
%==========================================================================


%==========================================================================
function plotFittingResult
if FitOpt.Plot && UserCommand~=99
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
%==========================================================================


%==========================================================================
function deleteSetButtonCallback(object,src,event)
global FitData
h = findobj('Tag','SetListBox');
idx = h.Value;
str = h.String;
nSets = numel(str);
if nSets>0
  ID = sscanf(str{idx},'%d');
  for k = numel(FitData.FitSets):-1:1
    if FitData.FitSets(k).ID==ID
      FitData.FitSets(k) = [];
    end
  end
  if idx>length(FitData.FitSets), idx = length(FitData.FitSets); end
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
return
%==========================================================================


%==========================================================================
function deleteSetListKeyPressFcn(object,event)
if strcmp(event.Key,'delete')
  deleteSetButtonCallback(object,gco,event);
  displayFitSet
end
return
%==========================================================================


%==========================================================================
function setListCallback(object,src,event)
  displayFitSet
return
%==========================================================================


%==========================================================================
function displayFitSet
global FitData
h = findobj('Tag','SetListBox');
idx = h.Value;
str = h.String;
if ~isempty(str)
  ID = sscanf(str{idx},'%d');
  k = find([FitData.FitSets.ID]==ID);  
  if k>0
    fitset = FitData.FitSets(k);
    values = fitset.bestvalues;
    
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

return
%==========================================================================


%==========================================================================
function exportSetButtonCallback(object,src,event)
global FitData
h = findobj('Tag','SetListBox');
v = h.Value;
s = h.String;
ID = sscanf(s{v},'%d');
idx = find([FitData.FitSets.ID]==ID);
varname = sprintf('fit%d',ID);
fitSet = rmfield(FitData.FitSets(idx),'bestx');
assignin('base',varname,fitSet);
fprintf('Fit set %d assigned to variable ''%s''.\n',ID,varname);
evalin('base',varname);
return
%==========================================================================


%==========================================================================
function selectAllButtonCallback(object,src,event)
h = findobj('Tag','ParameterTable');
d = h.Data;
d(:,1) = {true};
set(h,'Data',d);
return
%==========================================================================


%==========================================================================
function selectNoneButtonCallback(object,src,event)
h = findobj('Tag','ParameterTable');
d = h.Data;
d(:,1) = {false};
set(h,'Data',d);
return
%==========================================================================


%==========================================================================
function selectInvButtonCallback(object,src,event)
h = findobj('Tag','ParameterTable');
d = h.Data;
for k=1:size(d,1)
  d{k,1} = ~d{k,1};
end
set(h,'Data',d);
return
%==========================================================================


%==========================================================================
function sortIDSetButtonCallback(object,src,event)
global FitData
ID = [FitData.FitSets.ID];
[~,idx] = sort(ID);
FitData.FitSets = FitData.FitSets(idx);
refreshFitsetList(0);
return
%==========================================================================


%==========================================================================
function sortRMSDSetButtonCallback(object,src,event)
global FitData
rmsd = [FitData.FitSets.rmsd];
[~,idx] = sort(rmsd);
FitData.FitSets = FitData.FitSets(idx);
refreshFitsetList(0);
return
%==========================================================================


%==========================================================================
function refreshFitsetList(idx)
global FitData
h = findobj('Tag','SetListBox');
nSets = numel(FitData.FitSets);
for k=1:nSets
  s{k} = sprintf('%d. rmsd %g (%s)',...
    FitData.FitSets(k).ID,FitData.FitSets(k).rmsd,FitData.FitSets(k).Target);
end
if nSets==0, s = {}; end
set(h,'String',s);
if idx>0, set(h,'Value',idx); end
if idx==-1, set(h,'Value',numel(s)); end

if nSets>0, state = 'on'; else, state = 'off'; end
set(findobj('Tag','deleteSetButton'),'Enable',state);
set(findobj('Tag','exportSetButton'),'Enable',state);
set(findobj('Tag','sortIDSetButton'),'Enable',state);
set(findobj('Tag','sortRMSDSetButton'),'Enable',state);

displayFitSet;
return
%==========================================================================


%==========================================================================
function saveFitsetCallback(object,src,event)
global FitData
FitData.lastSetID = FitData.lastSetID+1;
FitData.currFitSet.ID = FitData.lastSetID;
if ~isfield(FitData,'FitSets') || isempty(FitData.FitSets)
  FitData.FitSets(1) = FitData.currFitSet;
else
  FitData.FitSets(end+1) = FitData.currFitSet;
end
refreshFitsetList(-1);
return
%==========================================================================


%==========================================================================
function tableEditCallback(hTable,callbackData)
global FitData

hTable = callbackData.Source;

% Get row and column index of edited table cell
ridx = callbackData.Indices(1);
cidx = callbackData.Indices(2);

% Return unless it's the last two columns
if cidx<5, return; end

% Revert if user-entered string does not cleanly convert to a scalar.
newval = str2num(callbackData.EditData);
if isempty(newval) || numel(newval)~=1
  hTable.Data{ridx,cidx} = callbackData.PreviousData;
  warning('Input is not a number.');
  return
end

% Get lower and upper bounds of interval from table
if cidx==5
  lower = newval;
  upper = str2num(hTable.Data{ridx,6});
elseif cidx==6
  lower = str2num(hTable.Data{ridx,5});
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
if FitData.nSystems>1
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

return
