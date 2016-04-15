% levelsplot    plot energy level diagram 
%
%  levelsplot(Sys,Ori,B);
%  levelsplot(Sys,Ori,B,mwFreq);
%  levelsplot(Sys,Ori,B,mwFreq,Par);
%
%    Sys        spin system structure
%    Ori        2-element vector [phi theta]
%               Euler angles for the magnetic field
%               alternatively, either 'x', 'y' or 'z'
%    B          field range [Bmin Bmax] in mT
%               alternatively, just Bmax in mT
%    mwFreq     spectrometer frequency in GHz
%    Par        other parameters
%        nPoints:     Number of points
%        ColorThreshold: Coloring threshold. All transitions with
%                        relative intensity below this will be gray.
%                        Example: 0.01
%        PlotThreshold:  All transitions below with relative intensity
%                        below this value will not be plotted.
%
%  If mwFreq is given, resonances are drawn. Red
%  lines indicate allowed transitions, gray lines
%  forbidden ones. If the lines are terminated with
%  dots, the relative transition intensity is larger
%  than 1%.
%
%  Example:
%    Sys = struct('S',7/2,'g',[2 2 2],'D',[1 1 -2]*5e3);
%    levelsplot(Sys,[0;pi/3],[0 6e3],95);

function levelsplot(Sys,varargin)

Para = struct;
switch (nargin)
case 5
  [Ori,B,mwFreq,Para] = deal(varargin{:});
case 4
  [Ori,B,mwFreq] = deal(varargin{:});
case 3
  [Ori,B] = deal(varargin{:});
  mwFreq = inf;
case 2
  error('Third input (magnetic field range [Bmin Bmax]) missing.');
case 0
  help(mfilename);
  return;
otherwise
  error('Wrong number of input arguments!');
end

if (nargout<0), error('Not enough output arguments.'); end
if (nargout>0), error('Too many output arguments.'); end

if isstruct(Ori)
  error('Second input argument (Ori) wrong: can''t be a structure.');
end
if isstruct(B)
  error('Third input argument (B) wrong: can''t be a structure.');
end

if ~isfield(Para,'nPoints')
  Para.nPoints = 201;
end
if ~isfield(Para,'PlotThreshold')
  Para.PlotThreshold = 0.001;
end
if ~isfield(Para,'ColorThreshold')
  Para.ColorThreshold = 0.01;
end

if ischar(Ori)
  switch Ori
    case 'x', Ori = [0, pi/2];
    case 'y', Ori = [pi/2, pi/2];
    case 'z', Ori = [0, 0];
    case 'xy', Ori = [pi/4 pi/2];
    case 'xz', Ori = [0 pi/4];
    case 'yz', Ori = [pi/2 pi/4];
    case 'xyz', Ori = [pi/4 acos(1/sqrt(3))];
    otherwise
      error('Unknown value ''%s'' for orientation (2nd input argument).',Ori);
  end
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
  error('Ori must be a two-element or three-element array [phi theta] or [phi theta chi].');
end

switch numel(B)
  case 1
    B = [0 B];
    Bvec = linspace(B(1),B(2),Para.nPoints);
  case 2
    Bvec = linspace(B(1),B(2),Para.nPoints);
  otherwise
    Bvec = B;
end

E = levels(Sys,[phi theta chi],Bvec);

if max(Bvec)>=2000
  Bscale = 1e-3;
else
  Bscale = 1;
end

plot(Bvec*Bscale,E/1e3,'b');

AllowedColor = [1 0 0];
ForbiddenColor = [1 1 1]*0.8;

if ~isfield(Sys,'S')
  nElectrons = 1;
  Sys.S = 1/2;
else
  nElectrons = numel(Sys.S);
end

if ~isfield(Sys,'g');
  Sys.g = ones(1,nElectrons)*2;
end

if isfinite(mwFreq)
  Opt = struct('Threshold',0,'Freq2Field',0);
  Exp = struct('mwFreq',mwFreq,'Range',B([1 end]));
  Exp.CrystalOrientation = [phi theta chi];
  [resonFields,tp,w,Transitions] = resfields(Sys,Exp,Opt);
  if ~isempty(resonFields)
    if (nElectrons==1)
      % one electron spin: normalize amplitudes
      %n = ang2vec(phi,theta);
      %g = norm(diag(Sys.g)*n);
      %c = g*bmagn/planck/1e9;
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
      if tp(iF)<Para.PlotThreshold, continue; end
      H = F + G*resonFields(iF);
      E = sort(eig(H))/1e3;
      h = line(resonFields(iF)*[1 1]*Bscale,E(Transitions(iF,:)));
      Color = tp(iF)*AllowedColor + (1-tp(iF))*ForbiddenColor;
      set(h,'Color',Color);
      if tp(iF)>Para.ColorThreshold, set(h,'Marker','.'); end
    end
    
  else
    % no resonance fields
    xl = xlim; yl = ylim;
    h = text(xl(1),yl(1),' no resonances in range!');
    set(h,'Color','r','VerticalAl','bottom');
  end
end

axis tight;
xlabel('magnetic field (mT)');
ylabel('energy (GHz)');
if isfinite(mwFreq)
  title(sprintf('%g GHz',mwFreq));
end

return