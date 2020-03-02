% levelsplot    plot energy level (Zeeman) diagram 
%
%  levelsplot(Sys,Ori,B)
%  levelsplot(Sys,Ori,B,mwFreq)
%  levelsplot(Sys,Ori,B,mwFreq,Par)
%
%    Sys        spin system structure
%    Ori        orientation of magnetic field in molecular frame
%               - string: 'x','y','z','xy','xz','yz', or 'xyz'
%               - 2-element vector [phi theta]
%               - 3-element vector [phi theta chi]
%    B          field range, in mT; either Bmax, [Bmin Bmax], or a full vector
%    mwFreq     spectrometer frequency, in GHz
%    Par        other parameters
%      Units           energy units for plotting, 'GHz' or 'cm^-1' or 'eV'
%      nPoints         Number of points
%      ColorThreshold  Coloring threshold. All transitions with relative
%                      intensity below this will be gray. Example: 0.05
%      PlotThreshold   All transitions below with relative intensity
%                      below this value will not be plotted. Example: 0.005
%
%  If mwFreq is given, resonances are drawn. Red lines indicate allowed
%  transitions, gray lines forbidden ones. If the lines are terminated with
%  dots, the relative transition intensity is larger than 1%.
%
%  Example:
%    Sys = struct('S',7/2,'g',2,'D',5e3);
%    levelsplot(Sys,'xy',6e3,95);

function levelsplot(Sys,varargin)

if nargin==0, help(mfilename); return; end 

% Default parameters
Ori = 'z';
B = [0 1400]; % mT
mwFreq = inf; % GHz
Par = struct;

% Parse input arguments
switch (nargin)
case 5
  Ori = varargin{1};
  B = varargin{2};
  mwFreq = varargin{3};
  Par = varargin{4};
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
if ~isnumeric(mwFreq) || numel(mwFreq)~=1 || ~isreal(mwFreq)
  error('Fourth input argument (mwFreq) must be a single number.');
end

% Supply option defaults
if ~isstruct(Par)
  error('Fifth input (parameters) must be a structure.');
end
if ~isfield(Par,'nPoints')
  Par.nPoints = 201;
end
if ~isfield(Par,'PlotThreshold')
  Par.PlotThreshold = 0.001;
end
if ~isfield(Par,'ColorThreshold')
  Par.ColorThreshold = 0.01;
end
Par.AllowedColor = [1 0 0];
Par.ForbiddenColor = [1 1 1]*0.8;
if ~isfield(Par,'Units')
  Par.Units = 'GHz';
end

% Convert string in Ori input to angles
if ischar(Ori)
  n = letter2vec(Ori);
  Ori = vec2ang(n);
end

if (numel(Ori)==2)
  phi = Ori(1);
  theta = Ori(2);
  chi = 0;
elseif (numel(Ori)==3)
  phi = Ori(1);
  theta = Ori(2);
  chi = Ori(3);
else
  error('Orientation (2nd input argument) must be a string (''x'',''y'',''z'',''xy'',''xz'',''yz'',''xyz''), a two-element array [phi theta], or a three-element array [phi theta chi].');
end

switch numel(B)
  case 1
    B = [0 B];
    Bvec = linspace(B(1),B(2),Par.nPoints);
  case 2
    Bvec = linspace(B(1),B(2),Par.nPoints);
  otherwise
    Bvec = B;
end

E_MHz = levels(Sys,[phi theta chi],Bvec);
E = unitconvert(E_MHz,Par.Units);

% Set horizontal and vertical units, scaling and labels
if max(Bvec)>=2000 % mT
  Bscale = 1e-3;  % use tesla for plotting
  xLabel = 'magnetic field (T)';
else
  Bscale = 1; % use millitesla for plotting
  xLabel = 'magnetic field (mT)';
end

switch Par.Units
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

h = plot(Bvec*Bscale,E*Escale,'b');
set(h,'Color',[0 0.4470 0.7410]);
axis tight
xlabel(xLabel);
ylabel(yLabel);

computeResonances = isfinite(mwFreq);

if computeResonances
  if ~isfield(Sys,'S')
    nElectrons = 1;
    Sys.S = 1/2;
  else
    nElectrons = numel(Sys.S);
  end
  
  if ~isfield(Sys,'g')
    Sys.g = ones(1,nElectrons)*2;
  end
  
  Opt = struct('Threshold',0,'Freq2Field',0);
  Exp = struct('mwFreq',mwFreq,'Range',B([1 end]));
  Exp.CrystalOrientation = [phi theta chi];
  [resonFields,tp,w,Transitions] = resfields(Sys,Exp,Opt); %#ok<ASGLU>
  tp = abs(tp); % to cover emissive transitions (in spin-polarized systems)
  if ~isempty(resonFields)
    if nElectrons==1
      % one electron spin: normalize amplitudes
      %n = ang2vec(phi,theta);
      %geff = norm(diag(Sys.g)*n);
      %c = geff*bmagn/planck/1e9;
      %tpMax = floor(Sys.S+0.5)*ceil(Sys.S+0.5)/4*c^2;
      %tp(tp>tpMax) = tpMax; % guard against round-off problems
      tpMax = max(tp);
    else
      tpMax = max(tp);
    end
    if (tpMax>0), tp = tp/tpMax; end
    
    % sort transitions according to intensity to ensure more intense
    % lines are plotted on top of less intense ones
    [tp,ix] = sort(tp);
    resonFields = resonFields(ix);
    Transitions = Transitions(ix,:);

    % compute and plot lower and upper energy levels of transitions
    n = ang2vec(phi,theta);
    [F,G] = sham(Sys,n);
    for iF = 1:numel(resonFields)
      if tp(iF)<Par.PlotThreshold, continue; end
      H = F + G*resonFields(iF);
      E_MHz = sort(eig(H));
      E = unitconvert(E_MHz,Par.Units);
      h = line(resonFields(iF)*[1 1]*Bscale,Escale*E(Transitions(iF,:)));
      Color = tp(iF)*Par.AllowedColor + (1-tp(iF))*Par.ForbiddenColor;
      set(h,'Color',Color);
      if tp(iF)>Par.ColorThreshold, set(h,'Marker','.'); end
      set(h,'ButtonDownFcn',sprintf('disp(''%d<->%d  transition probability = %0.4g'');',...
        Transitions(iF,1),Transitions(iF,2),tp(iF)));
    end
    
  else
    % no resonance fields
    xl = xlim; yl = ylim;
    h = text(xl(1),yl(1),' no resonances in field range');
    set(h,'Color','r','VerticalAl','bottom');
  end
end

if isfinite(mwFreq)
  switch Par.Units
    case 'GHz'
      str = sprintf('  %g GHz',mwFreq);
    case 'cm^-1'
      str = sprintf('  %0.3g cm^{-1} (%g GHz)',mwFreq*1e9/clight/100,mwFreq);
    case 'eV'
      if Escale==1e3
        str = sprintf('  %0.3g meV (%g GHz)',mwFreq*1e9*planck/evolt*1e3,mwFreq);
      else
        str = sprintf('  %0.3g eV (%g GHz)',mwFreq*1e9*planck/evolt,mwFreq);
      end
  end
  xl = xlim; yl = ylim;
  text(xl(1),yl(2),str,'VerticalAl','top');
end

return

function Eout = unitconvert(E_MHz,toUnit)

switch toUnit
  case 'GHz', Eout = E_MHz/1e3;
  case 'cm^-1', Eout = E_MHz*1e6/clight/100;
  case 'eV', Eout = planck*E_MHz*1e6/evolt;
  otherwise
    error('Unknown unit ''%s'' in Par.Units.',toUnit);
end
