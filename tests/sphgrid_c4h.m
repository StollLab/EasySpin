function ok = test()

% Test whether 1-octant closed-phi (D2h) and 1-octant open-phi (C4h) grids
% have the same number of triangles.

GridSize = [2 12 30];

for i = 1:numel(GridSize)
  [~,tri1] = sphgrid('D2h',GridSize(i));
  [~,tri2] = sphgrid('C4h',GridSize(i));
  ok(i) = numel(tri1.areas)==numel(tri2.areas);
end
