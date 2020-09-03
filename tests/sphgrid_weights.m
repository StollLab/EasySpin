function ok = test()

% Verify that weights always add up to 4*pi

nKnots = 10;

SymmGroups = {'C1','Ci','C2h','S6','C4h','C6h','D2h','Th',...
  'D3d','D4h','Oh','D6h','Dinfh','O3'};

for g = 1:numel(SymmGroups)
  grid = sphgrid(SymmGroups{g},nKnots);
  wsum(g) = sum(grid.weights);
end

ok = areequal(wsum,4*pi*ones(size(wsum)),1e-10,'rel');
