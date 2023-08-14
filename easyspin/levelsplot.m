% levelsplot    Plot energy level (Zeeman) diagram 
%
%  levelsplot(Sys,Ori,B)
%  levelsplot(Sys,Ori,B,mwFreq)
%  levelsplot(Sys,Ori,B,mwFreq,Opt)
%
%    Sys        spin system structure
%    Ori        orientation of magnetic field vector in molecular frame
%               - string: 'x','y','z','xy','xz','yz', or 'xyz'
%               - 2-element vector [phi theta] (radians)
%               orientation of lab frame in molecular frame
%               - 3-element vector [phi theta chi] (radians)
%    B          field range, in mT; either Bmax, [Bmin Bmax], or a full vector
%    mwFreq     spectrometer frequency, in GHz
%    Opt        options
%      Units           energy units for plotting, 'GHz' or 'cm^-1' or 'eV'
%      nPoints         Number of points
%      PlotThreshold   All transitions with relative intensity below
%                      this value will not be plotted. Example: 0.005
%      SlopeColor      true/false. Color energy level lines by their slope,
%                      (corresponding to their mS expectation value).
%                      Default is false.
%      StickSpectrum   true/false. Plot stick spectrum underneath Zeeman
%                      diagram. Default is false.
%
%  If mwFreq is given, resonances are drawn. Red lines indicate allowed
%  transitions, gray lines forbidden ones.
%
%  Example:
%    Sys = struct('S',7/2,'g',2,'D',5e3);
%    levelsplot(Sys,'xy',6e3,95);

function levelsplot(Sys,varargin)

if nargin==0, help(mfilename); return; end 

% Default parameters
Ori = 'z';
B = [0 1400];  % mT
mwFreq = inf;  % GHz
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
    fprintf('Third input (magnetic field range [Bmin Bmax]) missing; assuming [%d %d] mT.\n',B(1),B(2));
  case 1
    fprintf('Second input argument (orientation) missing; assuming ''%s'' (phi=0, theta=0, chi=0).\n',Ori);
    fprintf('Third input (magnetic field range [Bmin Bmax]) missing; assuming [%d %d] mT.\n',B(1),B(2));
  otherwise
    error('Too many input arguments!');
end

% Check number of output arguments
if nargout<0, error('Not enough output arguments.'); end
if nargout>0, error('Too many output arguments.'); end

% Check input parameters
if isstruct(Ori)
  error('Second input argument (Ori) wrong: can''t be a structure.');
end
if isstruct(B)
  error('Third input argument (B) wrong: can''t be a structure.');
end
if isempty(B) || numel(B)==1
  error('Input argument that specifies B range needs at least two elements, [Bmin Bmax].')
end
if ~isnumeric(mwFreq) || numel(mwFreq)~=1 || ~isreal(mwFreq)
  error('Fourth input argument (mwFreq) must be a single number.');
end
computeResonances = isfinite(mwFreq);

% Supply option defaults
if ~isstruct(Opt)
  error('Fifth input (options) must be a structure.');
end
if ~isfield(Opt,'nPoints')
  Opt.nPoints = 201;
end
if ~isfield(Opt,'PlotThreshold')
  Opt.PlotThreshold = 1e-6;
end
Opt.AllowedColor = [1 0 0];
Opt.ForbiddenColor = [1 1 1]*0.8;
if ~isfield(Opt,'Units')
  Opt.Units = 'GHz';
end
if ~isfield(Opt,'SlopeColor')
  Opt.SlopeColor = false;
end
if ~isfield(Opt,'StickSpectrum')
  Opt.StickSpectrum = false;
end

% Parse Ori (second input)
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

% Set ranges and units
%-------------------------------------------------------------------------------
switch numel(B)
  case 1
    B = [0 B];
    Bvec = linspace(B(1),B(2),Opt.nPoints);
  case 2
    Bvec = linspace(B(1),B(2),Opt.nPoints);
  otherwise
    Bvec = B;
end

% Calculate energy levels
E_MHz = levels(Sys,[phi theta chi],Bvec);
E = unit_convert(E_MHz,Opt.Units);
nLevels = size(E,2);

% Set horizontal and vertical units, scaling and labels
if max(Bvec)>=2000  % mT
  Bscale = 1e-3;  % use tesla for plotting
  xLabel = 'magnetic field (T)';
else
  Bscale = 1;  % use millitesla for plotting
  xLabel = 'magnetic field (mT)';
end

switch Opt.Units
  case 'GHz'
    yUnits = 'GHz';
    Escale = 1;
  case 'cm^-1'
    if max(abs(E(:)))<1e-2
      yUnits = '10^{-4} cm^{-1}';
      Escale = 1e4;
    else
      yUnits = 'cm^{-1}';
      Escale = 1;
    end
  case 'eV'
    if max(abs(E(:)))<1
      yUnits = 'meV';
      Escale = 1e3;
    else
      yUnits = 'eV';
      Escale = 1;
    end
  otherwise
    error('Unknown unit ''%s'' in Par.Units.',toUnit);
end
yLabel = ['energy (' yUnits ')'];

% Plot energy levels
%-------------------------------------------------------------------------------
if Opt.StickSpectrum && computeResonances
  subplot(4,1,[1 2 3]);
  cla
end
ax = gca;

if Opt.SlopeColor
  for iLevel = 1:nLevels
    col = abs(deriv(E(:,iLevel)));
    h = patch([Bvec(:); nan],[E(:,iLevel); nan],[col; nan],'EdgeColor','interp');
    h.Tag = 'level';
    h.UserData = iLevel;
  end
  colormap(flipud(parula));

else
  h = plot(Bvec*Bscale,E*Escale,'b');
  set(h,'Color',[0 0.4470 0.7410]);

  for iLevel = 1:nLevels
    h(iLevel).Tag = 'level';
    h(iLevel).UserData = iLevel;
  end
end
box on
axis tight
ylabel(yLabel);
set(gca,'Tag','diagram');
xl = xlim;
yl = ylim;
text(xl(1),yl(1),'','Tag','infotext','VerticalAlignment','bottom');

% Calculate and plot transitions
%-------------------------------------------------------------------------------
if computeResonances
  
  resfieldsOpt = struct('Threshold',0,'Freq2Field',0);
  Exp = struct('mwFreq',mwFreq,'Range',B([1 end]));
  Exp.SampleFrame = [-chi -theta -phi];
  [resonFields,intensity,~,Transitions] = resfields(Sys,Exp,resfieldsOpt);

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
    for iF = 1:numel(resonFields)
      if absintensity(iF)<Opt.PlotThreshold, continue; end

      H = H0 - muzL*resonFields(iF);
      E_MHz = sort(eig(H));
      E = unit_convert(E_MHz,Opt.Units);

      h = line(resonFields(iF)*[1 1]*Bscale,Escale*E(Transitions(iF,:)),'Tag','transition','LineWidth',1);
      h.UserData = [Transitions(iF,:) resonFields(iF) absintensity(iF)];
      transitionColor = absintensity(iF)*Opt.AllowedColor + (1-absintensity(iF))*Opt.ForbiddenColor;
      h.Color = transitionColor;
      h.ButtonDownFcn = @(src,~)fprintf('transition %d-%d:  relative intensity = %0.4g\n',...
        Transitions(iF,1),Transitions(iF,2),absintensity(iF));
    end

    if Opt.StickSpectrum
      xl = xlim;
      subplot(4,1,4);
      cla
      for iF = 1:numel(resonFields)
        hLine = line([1 1]*resonFields(iF)*Bscale,[0 intensity(iF)],'LineWidth',2,'Tag','line');
        hLine.UserData = [Transitions(iF,:) resonFields(iF) absintensity(iF)];
      end
      ylim([0 1.1]);
      xlim(xl);
      box on
    end
    
  else
    % no resonance fields
    xl = xlim;
    yl = ylim;
    h = text(xl(1),yl(1),' no resonances in field range');
    set(h,'Color','r','VerticalAl','bottom');
  end
end
xlabel(xLabel);

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
  str = sprintf('  %s (φ,θ,χ) = (%0.1f, %0.1f, %0.1f)°\n%s',...
                oristr,phi*180/pi,theta*180/pi,chi*180/pi,mwstr);
else
  str = sprintf('  %s (φ,θ) = (%0.1f, %0.1f)°\n%s',...
                oristr,phi*180/pi,theta*180/pi,mwstr);
end

xl = xlim(ax);
yl = ylim(ax);
text(ax,xl(1),yl(2),str,'VerticalAl','top');

% Activate mouseovers only if there is one levelsplot axes in the figure.
hAx = findobj(gcf,'Tag','diagram');
if numel(hAx)==1
  set(gcf,'WindowButtonMotionFcn',@windowButtonMotionFcn);
else
%  set(gcf,'WindowButtonMotionFcn','');
end

end

%-------------------------------------------------------------------------------
function Eout = unit_convert(E_MHz,toUnit)

switch toUnit
  case 'GHz', Eout = E_MHz/1e3;
  case 'cm^-1', Eout = E_MHz*1e6/clight/100;
  case 'eV', Eout = planck*E_MHz*1e6/evolt;
  otherwise
    error('Unknown unit ''%s'' in Par.Units.',toUnit);
end

end

%-------------------------------------------------------------------------------
function windowButtonMotionFcn(~,~,~)

% Obtain handle of object under mouse pointer
hObj = hittest();

% Construct information string for level, transition or spectral line
switch hObj.Tag
  case 'level'
    str = sprintf('level %d',hObj.UserData);
  case 'transition'
    str = sprintf('transition %d-%d: %0.2f mT, relative intensity %0.4f ',...
      hObj.UserData(1),hObj.UserData(2),hObj.UserData(3),hObj.UserData(4));
  case 'line'
    str = sprintf('transition %d-%d: %0.2f mT, relative intensity %0.4f ',...
      hObj.UserData(1),hObj.UserData(2),hObj.UserData(3),hObj.UserData(4));
  otherwise
    % Remove information if not over an object
    str = '';
end

% Display information
hParent = hObj.Parent;  % this is either the axes or the figure
hText = findobj(hParent,'Tag','infotext');
set(hText,'String',[' ' str]);

end
