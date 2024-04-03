% nucfrq2d  HYSCORE powder spectrum outline 
%
%   nucfrq2d(Sys,B0)
%   nucfrq2d(Sys,B0,tau)
%
%   Computes powder HYSCORE peaks (without amplitudes)
%   from spin system Sys at external magnetic field magnitude B0
%   (in mT) and plots the result.
%
%   alpha-beta correlations are in blue, beta-alpha correlations are in red.
%
%   Only S=1/2 systems are supported.
%
%   tau (in µs) specifies a vector of tau values.
%
%   If tau is given, a colored background indicating peak
%   suppression regions (blind spots) is shown. White indicates
%   no suppression, the darker the gray the stronger the suppression.
%
%   Example:
%
%      Sys = struct('Nucs','1H','g',[2 2 2]);
%      Sys.A = 3 + [-1 -1 2]*5;
%      nucfrq2d(Sys,350,0.120);   % field 350 mT, tau 120 ns

function nucfrq2d(sys,B0,tau)

if nargin==0, help(mfilename); return; end
if nargin<2 || nargin>3, error('Wrong number of input arguments!'); end
if nargout<0, error('Not enough output arguments.'); end
if nargout>0, error('Too many output arguments.'); end

% Check limitations.
if isfield(sys,'S')
  if sys.S>1/2
    error('Systems with S>1/2 are not supported.');
  end
end

% Internal options, mostly for plotting
%===============================================================================
Options.GridSize = 31;
Options.nPoints = 100;
Options.expand = 1.1;
Options.baColor = [1 0 0]*0.6;
Options.abColor = [0 0 1]*0.6;
Options.QuadraticAxes = 0;
%===============================================================================

if nargin<3, tau = 0; end

computeBlindSpots = all(tau>0) & ~Options.QuadraticAxes;

% Construct spin Hamiltonian and get state space dimension.
[H0,mux,muy,muz] = ham(sys);
N = size(H0,1);

% Construct masks for alpha and beta manifold transitions.
nn = N/2;
triup = triu(ones(nn),1);
BetaMask = logical(blkdiag(triup,zeros(nn)));
AlphaMask = logical(blkdiag(zeros(nn),triup));
nfreq = nnz(BetaMask);

% Set up orientational grid and triangulation.
thisSymm = hamsymm(sys);
[grid,tri] = sphgrid(thisSymm,Options.GridSize);
x = grid.vecs;
nOct = grid.nOctants;
nOri = size(x,2);

% Initialize to zero.
[FreqsAlpha,FreqsBeta] = deal(zeros(nOri,nfreq));

B = B0*x;

% Loop over all orientations and compute nuclear transitions.
for k = 1:nOri
  H = H0 - (B(1,k)*mux + B(2,k)*muy + B(3,k)*muz);
  E = sort(eig(H));
  EE = E(:,ones(1,N));
  EE = EE.' - EE;
  FreqsAlpha(k,:) = EE(AlphaMask).';
  FreqsBeta(k,:) = EE(BetaMask).';
end

if Options.QuadraticAxes
  FreqsAlpha = FreqsAlpha.^2;
  FreqsBeta = FreqsBeta.^2;
end

% Get maximum transition frequency, for plotting.
maxFrq = max(max(FreqsAlpha(:)),max(FreqsBeta(:)));

% Compute blind spot amplitude modulation.
if computeBlindSpots
  nPoints = 2*Options.nPoints+1;
  freq = Options.expand*maxFrq*linspace(-1,1,nPoints);
  Modulation = zeros(nPoints,nPoints);
  for k = 1:numel(tau)
    Modulation1D = abs(sin(pi*tau(k)*freq));
    Modulation = Modulation + Modulation1D.'*Modulation1D;
  end
  Modulation = Modulation/max(Modulation(:));
end

% Plotting
%===============================================================================

clf

% Blindspot pattern
if computeBlindSpots
  if Options.QuadraticAxes
    xy = freq.^2;
  else
    xy = freq;
  end
  pcolor(xy,xy,Modulation);
  ColMap = gray(256);
  colormap(ColMap(200:end,:))
  shading interp
end

hold on

% Diagonals and axes
line([0 0; 1 -1]*Options.expand*maxFrq,[0 0; 1 1]*Options.expand*maxFrq,'Color','k');
line([0 0],[0 1]*Options.expand*maxFrq,'Color','k');

% Antidiagonals at Larmor frequencies
LarmorFreqs = nmagn*B0*nucgval(sys.Nucs)/planck/1e9;
for k = 1:numel(LarmorFreqs)
  line([-2 0 2]*LarmorFreqs(k),[0 2 0]*LarmorFreqs(k),'Color','k','LineStyle','--');
end

baColor = Options.baColor;
abColor = Options.abColor;
if nOct<1
  % Axial case
  for i1 = 1:nfreq
    for i2 = 1:nfreq
      line(FreqsBeta(:,i1),FreqsAlpha(:,i2),'Color',baColor,'LineWidth',2);
      line(-FreqsBeta(:,i1),FreqsAlpha(:,i2),'Color',baColor,'LineWidth',2);
    end
  end
  for i1 = 1:nfreq
    for i2 = 1:nfreq
      line(FreqsAlpha(:,i1),FreqsBeta(:,i2),'Color',abColor,'LineWidth',2);
      line(-FreqsAlpha(:,i1),FreqsBeta(:,i2),'Color',abColor,'LineWidth',2);
    end
  end
else
  % Nonaxial case
  z = zeros(1,nOri);
  for i1 = 1:nfreq
    for i2 = 1:nfreq
      trisurf(tri.idx,FreqsBeta(:,i1),FreqsAlpha(:,i2),z,'EdgeColor',baColor,'FaceColor',baColor);
      trisurf(tri.idx,-FreqsBeta(:,i1),FreqsAlpha(:,i2),z,'EdgeColor',baColor,'FaceColor',baColor);
    end
  end
  for i1 = 1:nfreq
    for i2 = 1:nfreq
      trisurf(tri.idx,FreqsAlpha(:,i1),FreqsBeta(:,i2),z,'EdgeColor',abColor,'FaceColor',abColor);
      trisurf(tri.idx,-FreqsAlpha(:,i1),FreqsBeta(:,i2),z,'EdgeColor',abColor,'FaceColor',abColor);
    end
  end
end

axis equal
box on
axis(Options.expand*maxFrq*[-1 1 0 1]);
hold off
if Options.QuadraticAxes
  xlabel('\nu_1^2 (MHz^2)');
  ylabel('\nu_2^2 (MHz^2)');
else
  xlabel('\nu_1 (MHz)');
  ylabel('\nu_2 (MHz)');
end
set(gca,'Layer','top');

end
