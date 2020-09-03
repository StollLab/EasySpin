function ok = test()

% Syntax check

SymmGroup = {'D2h','C2h','Ci','Dinfh'};
nKnots = 10;
ok = true;
for g = 1:numel(SymmGroup)
  grid = sphgrid(SymmGroup{g},nKnots);
  ok(g) = isfield(grid,'phi') && ...
          isfield(grid,'theta') && ...
          isfield(grid,'vecs') && ...
          isfield(grid,'weights');
end

ok = all(ok);
