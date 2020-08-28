% sphtri  Triangulation of standard EasySpin spherical grid
%
%   tri = sphtri(Symmetry,n)
%
%   Computes the triangulation of the standard EasySpin spherical grid.
%
%   Inputs:
%     n        .. number of knots along a quarter of a meridian
%     Symmetry .. point group
%
%   Outputs:
%     tri      .. Nx3 matrix containing indices into the grid vector.
%
%   The triangulation of point groups with open phi intervals (C4h,C6h,S6,C2h)
%   is done after closing the phi interval.

% Undocumented:
%   tri = sphtri(nOctants,n)
%
%   corresponds to
%
%   nOctants = 1: Symmetry = 'D2h';
%   nOctants = 2: Symmetry = 'C2h';
%   nOctants = 4: Symmetry = 'Ci';
%   nOctants = 8: Symmetry = 'C1';

function Triangulation = sphtri(Symmetry,nKnots,Options)

if nargin==0, help(mfilename); return; end

if nargin<3, Options = ''; end

if nKnots<2
  error('Cannot compute triangulation with n<2.');
end

DebugMode = false;

if ~ischar(Options)
  error('Third input argument must be a string.');
end

explicitClosedPhi = any(Options=='f');

switch (Symmetry)
  case 1, Symmetry = 'D2h';
  case 2, Symmetry = 'C2h';
  case 4, Symmetry = 'Ci';
  case 8, Symmetry = 'C1';
end

switch (Symmetry)
  case {'D6h','D4h','Oh','D3d','Th','D2h',...
      'C4h','C6h'}, SymmID = 1;
  case 'Ci', SymmID = 2;
  case 'C1', SymmID = 3;
  case {'C2h','S6'}, SymmID = 4;
  case 'Dinfh', SymmID = 5;
  otherwise
    error('Symmetry not supported!');
end

% Cache previously computed triangulations.
persistent triknown
if ~isempty(triknown)
  idx = find((triknown.nKnots==nKnots) & ...
             (triknown.SymmID==SymmID) & ...
             (triknown.explicitClosedPhi==explicitClosedPhi));
  if ~isempty(idx)
    Triangulation = triknown.Triangulation{idx};
    if DebugMode
      disp(sprintf('sphtri: retrieve, nKnots %d, SymmID %d, expl %d.',nKnots,SymmID));
    end
    return
  end
else
  if DebugMode
    disp('sphtri: init.');
  end
  triknown.SymmID = [];
  triknown.nKnots = [];
  triknown.explicitClosedPhi = [];
  triknown.Triangulation = {};
end
if DebugMode
  disp('sphtri: compute.');
end

switch SymmID

  case 1 %'D6h','D4h','Oh','D3d','Th','D2h',... % closed phi interval
      %'C4h','C6h'} % open phi interval
    % coding idea of David Goodmanson, comp.soft-sys.matlab
    a = 1:nKnots*(nKnots+1)/2;
    a((2:nKnots+1).*(1:nKnots)/2) = [];
    k = 1:(nKnots-1)*(nKnots-2)/2;
    b = a(k);
    Triangulation = [[1:nKnots*(nKnots-1)/2; a; a+1],[b; 2*(b+1)-k; b+1]].';

  case 2 %'Ci' % 4 octants
    % Same number of triangles 4*(nKnots-1)^2 for both open and closed grid!
    if explicitClosedPhi
      [phx,thx] = sphgrid('Ci',nKnots,'f');
      phx = phx(:); thx = thx(:);
      q = 1/2;
      phx = q*pi+(phx-q*pi).*(1-0.01*thx);
    else
      [phx,thx] = sphgrid('Ci',nKnots);
      phx = phx(:); thx = thx(:);
      phx = phx + rand(size(phx))*(0.001*pi/180); % remove collinearity
      q = 1;
    end
    x = thx.*cos(q*phx);
    y = thx.*sin(q*phx);
    Triangulation = delaunayn([x,y]);
    % Same number of triangles 4*(nKnots-1)^2 for both
    % open and closed grid!

  case 3 %'C1'
    [phx,thx] = sphgrid('Ci',nKnots);
    phx = phx + rand(size(phx))*(0.001*pi/180); % remove collinearity
    Triangulation = delaunay(thx.*cos(phx), thx.*sin(phx));
    Triangulation1 = (4*nKnots^2 - 8*nKnots + 6) + 1 - Triangulation;
    idx2 = Triangulation > length(thx)-4*(nKnots-1);
    Triangulation1(idx2) = Triangulation(idx2);
    Triangulation = [Triangulation; Triangulation1];

  case 4 % {'C2h','S6'} % 2 octants, periodic
    [phx,thx] = sphgrid(Symmetry,nKnots,'f');
    phx = phx(:); thx = thx(:);
    phx2 = max(phx)/2;
    phx = phx2 + (phx-phx2).*(1-thx/500); % make fully convex
    Triangulation = delaunayn([thx.*cos(phx), thx.*sin(phx)]);
    %idx = (1:nKnots).^2;
    %for k=2:nKnots
    %  Triangulation(Triangulation==idx(k)) = idx(k-1)+1;
    %end

  case 5 %'Dinfh'
    Triangulation = [];
end

Triangulation = sort(Triangulation,2);
Triangulation = uint32(Triangulation);

n = numel(triknown.SymmID)+1;
triknown.Triangulation{n} = Triangulation;
triknown.nKnots(n) = nKnots;
triknown.SymmID(n) = SymmID;
triknown.explicitClosedPhi(n) = explicitClosedPhi;

return
