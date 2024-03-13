% levelsplot    Plot energy level (Zeeman) diagram 
%
%  levelsplot(Sys,Ori,B)
%  levelsplot(Sys,Ori,B,mwFreq)
%  levelsplot(Sys,Ori,B,mwFreq,Opt)
%
%    Sys        spin system structure
%    Ori        (a) orientation of magnetic field vector in molecular frame
%               - string: 'x','y','z','xy','xz','yz', or 'xyz'
%               - 2-element vector [phi theta] (radians)
%               (b) orientation of lab frame in molecular frame
%               - 3-element vector [phi theta chi] (radians)
%    B          field range, in mT; either Bmax, [Bmin Bmax], or a full vector
%    mwFreq     spectrometer frequency, in GHz
%    Opt        options
%      Units           energy units for plotting, 'GHz' or 'cm^-1' or 'eV'
%      nPoints         number of points
%      PlotThreshold   all transitions with relative intensity below
%                      this value will not be plotted. Example: 0.005
%      SlopeColor      true/false (default false). Color energy level lines
%                      by their slope, (corresponding to their mS expectation
%                      value).
%      StickSpectrum   true/false (default false). Plot stick spectrum underneath
%                      Zeeman diagram. Default is false.
%
%  If mwFreq is given, resonances are drawn. Red lines indicate allowed
%  transitions, gray lines forbidden ones. Hovering with the cursor over
%  the lines displays intensity information.
%
%  Example:
%    Sys = struct('S',7/2,'g',2,'D',5000);
%    levelsplot(Sys,'xy',[0 6000],95);

function levelsplot(Sys,varargin)

if nargin==0, help(mfilename); return; end 

% Default values for input arguments
Ori = 'z';
B = [0 1400];  % mT
mwFreq = [];  % GHz
Opt = struct;

% Parse input arguments
%-------------------------------------------------------------------------------
switch nargin
  case 5
    Ori = varargin{1};
    B = varargin{2};
    mwFreq = varargin{3};
    Opt = varargin{4};
  case 4
    Ori = varargin{1};
    B = varargin{2};
    mwFreq = varargin{3};
  case 3
    Ori = varargin{1};
    B = varargin{2};
  case 2
    Ori = varargin{1};
    fprintf('Third input argument (magnetic field range [Bmin Bmax]) missing; assuming [%d %d] mT.\n',B(1),B(2));
  case 1
    fprintf('Second input argument (orientation) missing; assuming ''%s''.\n',Ori);
    fprintf('Third input argument (magnetic field range [Bmin Bmax]) missing; assuming [%d %d] mT.\n',B(1),B(2));
  otherwise
    error('Too many input arguments! At most 5 (Sys,Ori,B,mwFreq,Opt) are possible.');
end

% Check number of output arguments
%-------------------------------------------------------------------------------
if nargout~=0, error('levelsplot does not return output arguments.'); end

% Parse Options (fifth input argument)
%-------------------------------------------------------------------------------
if ~isstruct(Opt)
  error('Fifth input argument (options) must be a structure.');
end
if ~isfield(Opt,'Units')
  Opt.Units = 'GHz';
end
if ~isfield(Opt,'nPoints')
  Opt.nPoints = 201;
end
if ~isfield(Opt,'PlotThreshold')
  Opt.PlotThreshold = 1e-6;
end
Opt.AllowedColor = [1 0 0];
Opt.ForbiddenColor = [1 1 1]*0.8;
if ~isfield(Opt,'SlopeColor')
  Opt.SlopeColor = false;
end
if ~isfield(Opt,'StickSpectrum')
  Opt.StickSpectrum = false;
end

% Parse Ori (second input argument)
%-------------------------------------------------------------------------------
if isstruct(Ori)
  error('Second input argument (Ori) wrong: can''t be a structure.');
end
if ischar(Ori)
  n = letter2vec(Ori);
  [phi,theta] = vec2ang(n);
  chi = 0;
elseif numel(Ori)==2
  phi = Ori(1);
  theta = Ori(2);
  chi = 0;
elseif numel(Ori)==3
  phi = Ori(1);
  theta = Ori(2);
  chi = Ori(3);
else
  error('Orientation (2nd input argument) must be a string (''x'',''y'',''z'',''xy'',''xz'',''yz'',''xyz''), a two-element array [phi theta], or a three-element array [phi theta chi].');
end

% Parse B (third input argument)
%-------------------------------------------------------------------------------
if isstruct(B)
  error('Third input argument (B) wrong: can''t be a structure.');
end
if isempty(B) || numel(B)==1
  error('Input argument that specifies B range needs at least two elements, [Bmin Bmax].')
end
switch numel(B)
  case 2
    Bvec = linspace(B(1),B(2),Opt.nPoints);
  otherwise
    Bvec = B;
end
% Set horizontal and vertical units, scaling and labels
if max(Bvec)>=2000  % mT
  Bscale = 1e-3;  % use tesla for plotting
  fieldUnit = 'T';
else
  Bscale = 1;  % use millitesla for plotting
  fieldUnit = 'mT';
end

% Parse mwFreq (fourth input argument)
%-------------------------------------------------------------------------------
if ~isempty(mwFreq)
  if ~isnumeric(mwFreq) || numel(mwFreq)~=1 || ~isreal(mwFreq)
    error('Fourth input argument (mwFreq) must be a single number.');
  end
end
computeResonances = ~isempty(mwFreq);
if ~computeResonances && Opt.StickSpectrum
  warning('Cannot plot stick spectrum, since microwave frequency is missing.');
end

% Calculate energy levels
%-------------------------------------------------------------------------------
E_MHz = levels(Sys,[phi theta chi],Bvec);
E = unit_convert(E_MHz,Opt.Units);
nLevels = size(E,2);

switch Opt.Units
  case 'GHz'
    energyUnit = 'GHz';
    Escale = 1;
  case 'cm^-1'
    if max(abs(E(:)))<1e-2
      energyUnit = '10^{-4} cm^{-1}';
      Escale = 1e4;
    else
      energyUnit = 'cm^{-1}';
      Escale = 1;
    end
  case 'eV'
    if max(abs(E(:)))<1
      energyUnit = 'meV';
      Escale = 1e3;
    else
      energyUnit = 'eV';
      Escale = 1;
    end
  otherwise
    % This was already caught by unit_convert()
    error('Unsupported unit ''%s'' in Opt.Units.',toUnit);
end

% Plot energy levels
%-------------------------------------------------------------------------------
if Opt.StickSpectrum && computeResonances
  subplot(4,1,[1 2 3]);
end
cla
hLevelsAxes = gca;
linecolor = [0 0.4470 0.7410];

if Opt.SlopeColor
  for iLevel = nLevels:-1:1
    col = abs(deriv(E(:,iLevel)));
    hLevels(iLevel) = patch([Bvec(:)*Bscale; NaN],[E(:,iLevel); NaN],[col; NaN],'EdgeColor','interp');
  end
  colormap(flipud(parula));
else
  hLevels = plot(Bvec*Bscale,E*Escale,'b');
  set(hLevels,'Color',linecolor);
end
for iLevel = 1:nLevels
  hLevels(iLevel).Tag = 'level';
  hLevels(iLevel).UserData = iLevel;
end

try
  % Adjust axes interactivity and axes toolbar
  hLevelsAxes.Interactions = [];
  axtoolbar(hLevelsAxes,{'export','pan','zoomin','zoomout','restoreview'});
catch
  % pre-R2018b
end

box on
axis tight
ylabel(sprintf('energy (%s)',energyUnit));
set(gca,'Tag','diagram');
xl = xlim;
yl = ylim;
text(xl(1),yl(1),'','Tag','infotext','VerticalAlignment','bottom');

% Calculate and plot transitions if requested
%-------------------------------------------------------------------------------
if computeResonances
  
  % Calculate resonance fields
  resfieldsOpt = struct('Threshold',0,'Freq2Field',0);
  Exp = struct('mwFreq',mwFreq,'Range',B([1 end]));
  Exp.SampleFrame = [-chi -theta -phi];
  [resonFields,intensity,~,Transitions] = resfields(Sys,Exp,resfieldsOpt);

  % Plot transitions at resonance fields
  if ~isempty(resonFields)
    tpMax = max(abs(intensity));
    if tpMax>0, intensity = intensity/tpMax; end
    absintensity = abs(intensity);
    % sort transitions according to intensity to ensure more intense
    % lines are plotted on top of less intense ones
    [absintensity,idx] = sort(absintensity);
    resonFields = resonFields(idx);
    Transitions = Transitions(idx,:);
    intensity = intensity(idx);

    % compute and plot lower and upper energy levels of transitions
    zL = ang2vec(phi,theta);  % lab z direction in molecular frame representation
    [H0,muzL] = ham(Sys,zL);
    for iF = numel(resonFields):-1:1
      if absintensity(iF)<Opt.PlotThreshold, continue; end

      H = H0 - muzL*resonFields(iF);
      E_MHz = sort(eig(H));
      E = unit_convert(E_MHz,Opt.Units);

      h = line(resonFields(iF)*[1 1]*Bscale,Escale*E(Transitions(iF,:)),'Tag','transition','LineWidth',1);
      h.UserData = [Transitions(iF,:) resonFields(iF) absintensity(iF)];
      transitionColor = absintensity(iF)*Opt.AllowedColor + (1-absintensity(iF))*Opt.ForbiddenColor;
      h.Color = transitionColor;
      h.ButtonDownFcn = @(src,~)fprintf('transition %d-%d:  %g mT, relative intensity = %0.4g\n',...
        Transitions(iF,1),Transitions(iF,2),resonFields(iF),absintensity(iF));
    end

    if Opt.StickSpectrum
      xl = xlim;
      subplot(4,1,4);
      cla
      yline(0);
      for iF = 1:numel(resonFields)
        hLine = line([1 1]*resonFields(iF)*Bscale,[0 intensity(iF)],'Color',linecolor,'LineWidth',2,'Tag','line');
        hLine.UserData = [Transitions(iF,:) resonFields(iF) absintensity(iF)];
      end
      if any(intensity<0)
        ylim([-1 1]*1.1);
      else
        ylim([0 1.1]);
      end
      xlim(xl);
      box on
      set(gca,'FontSize',hLevelsAxes.FontSize)
      hStickSpectrum = gca;
      hStickSpectrum.UserData = linkprop([hStickSpectrum hLevelsAxes],'XLim');
    end
    
  else
    % no resonance fields
    xl = xlim;
    yl = ylim;
    h = text(xl(1),yl(1),' no resonances in field range');
    set(h,'Color','r','VerticalAl','bottom');
  end
end
xlabel(sprintf('magnetic field (%s)',fieldUnit));

% Display orientation and microwave frequency (if given)
%-------------------------------------------------------------------------------
if isfinite(mwFreq)
  switch Opt.Units
    case 'GHz'
      mwstr = sprintf('  %g GHz\n',mwFreq);
    case 'cm^-1'
      mwstr = sprintf('  %0.3g cm^{-1} (%g GHz)\n',mwFreq*1e9/clight/100,mwFreq);
    case 'eV'
      if Escale==1e3
        mwstr = sprintf('  %0.3g meV (%g GHz)\n',mwFreq*1e9*planck/evolt*1e3,mwFreq);
      else
        mwstr = sprintf('  %0.3g eV (%g GHz)\n',mwFreq*1e9*planck/evolt,mwFreq);
      end
  end
else
  mwstr = '';
end

if ischar(Ori)
  oristr = ['''' Ori '''' ', '];
  chiGiven = false;
else
  oristr = '';
  chiGiven = numel(Ori)==3;
end
if chiGiven
  oristr = sprintf('  %s (φ,θ,χ) = (%0.1f, %0.1f, %0.1f)°\n%s',...
                oristr,phi*180/pi,theta*180/pi,chi*180/pi,mwstr);
else
  oristr = sprintf('  %s (φ,θ) = (%0.1f, %0.1f)°\n%s',...
                oristr,phi*180/pi,theta*180/pi,mwstr);
end

xl = xlim(hLevelsAxes);
yl = ylim(hLevelsAxes);
text(hLevelsAxes,xl(1),yl(2),oristr,'VerticalAl','top');

% Activate mouseovers (callback function handles multiple axes)
%-------------------------------------------------------------------------------
set(gcf,'WindowButtonMotionFcn',@windowButtonMotionFcn);

end
%===============================================================================


%-------------------------------------------------------------------------------
function Eout = unit_convert(E_MHz,toUnit)

switch toUnit
  case 'GHz', Eout = E_MHz/1e3;
  case 'cm^-1', Eout = unitconvert(E_MHz,'MHz->cm^-1');
  case 'eV', Eout = unitconvert(E_MHz,'MHz->eV');
  otherwise
    error('Unsupported unit ''%s'' in Opt.Units.',toUnit);
end

end

%-------------------------------------------------------------------------------
function windowButtonMotionFcn(~,~,~)

persistent hPrevLine

hoverLineWidth = 2;
defaultLineWidth = 0.5;

% Obtain handle of object under mouse pointer
hObj = hittest(); % hittest() is an undocumented built-in MATLAB function :-/
hParent = hObj.Parent;  % this is either the axes or the figure

% Remove any previous thick lines
if ~isempty(hPrevLine) && ~strcmp(hPrevLine.Tag,'line')
  set(hPrevLine,'LineWidth',defaultLineWidth);
  hPrevLine = [];
end

% Construct information string for level, transition or spectral line,
% using UserData of object under mouse pointer
switch hObj.Tag
  case 'level'
    infostr = sprintf(' level %d',hObj.UserData);
    hObj.LineWidth = hoverLineWidth;
    hPrevLine = hObj;
  case 'transition'
    infostr = sprintf(' transition %d-%d: %0.2f mT, relative intensity %0.4f ',...
      hObj.UserData(1),hObj.UserData(2),hObj.UserData(3),hObj.UserData(4));
    hObj.LineWidth = hoverLineWidth;
    hPrevLine = hObj;
  case 'line'
    infostr = sprintf(' transition %d-%d: %0.2f mT, relative intensity %0.4f ',...
      hObj.UserData(1),hObj.UserData(2),hObj.UserData(3),hObj.UserData(4));
    hObj.LineWidth = hoverLineWidth;
    hPrevLine = hObj;
  otherwise
    % Remove information if not over a relevant object
    infostr = '';
end

% Display information
hText = findobj(hParent,'Tag','infotext');
set(hText,'String',infostr);

% Update position (for infotext to remain visible after zooming in)
if ~isempty(hText) && ~isempty(infostr)
  hText.Position(1:2) = [hParent.XLim(1) hParent.YLim(1)];
end

end
