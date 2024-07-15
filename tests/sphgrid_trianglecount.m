function ok = test()

% Assert that the correct number of triangles is returned for every point
% group symmetry.

% Test all point group symmetries
Symmetry = {'O3','Dinfh',...
  'D6h','D4h','Oh','D3d','Th','D2h','C4h','C6h',...
  'C2h','S6',...
  'Ci',...
  'C1'};

% Test several grid sizes
gridsize = [3 12 30];

for g = 1:numel(gridsize)
  nT1 = (gridsize(g)-1)^2;  % number of triangles in one octant
  nTriangles_ref = nT1*[0 0 1 1 1 1 1 1 1 1 2 2 4 8];
  for s = numel(Symmetry):-1:1
    [~,t] = sphgrid(Symmetry{s},gridsize(g));
    nTriangles(s) = numel(t.areas);
  end
  ok(g) = all(nTriangles==nTriangles_ref);
end
