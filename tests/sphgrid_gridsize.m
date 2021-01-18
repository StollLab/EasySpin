function ok = test(opt)

% Explicitly check grid size for all point groups

gridsize = containers.Map;
gridsize('C1') = @(m) m*(m-1)/2*4 + 1 + (m-1)*(m-2)/2*4 + 1; % 8 octants, open
gridsize('Ci') = @(m) m*(m-1)/2*4 + 1; % 4 octants, open

f = @(m) m*(m-1)/2*2 + 1;
gridsize('C2h') = f; % 2 octants, open
gridsize('S6') = f; % 2 octants, open

f = @(m) m*(m-1)/2 + 1;
gridsize('C4h') = f; % 1 octant, open
gridsize('C6h') = f; % 1 octant, open

f = @(m) (m+1)*m/2;
gridsize('D2h') = f; % 1 octant, closed
gridsize('D4h') = f; % 1 octant, closed
gridsize('D6h') = f; % 1 octant, closed
gridsize('Oh') = f; % 1 octant, closed
gridsize('Th') = f; % 1 octant, closed
gridsize('D3d') = f; % 1 octant, closed

GridSize = 12;

pg = keys(gridsize);
for idx = 1:length(gridsize)
  pg_ = pg{idx};
  grid = sphgrid(pg_,GridSize);
  n = numel(grid.weights);
  fun = gridsize(pg_);
  n0 = fun(GridSize);
  ok(idx) = n==n0;
  if opt.Display
    fprintf('%5s:  %d  %d\n',pg_,n,n0);
  end
end
