function ok = test()

% Explicitly check very small grid over full sphere

phi = [0 0 1/2 1 3/2 0]*pi;
theta = [0 1/2 1/2 1/2 1/2 1]*pi;
vecs = [0 0 1; 1 0 0; 0 1 0; -1 0 0; 0 -1 0; 0 0 -1].';
weights = [1.8403 2.2214 2.2214 2.2214 2.2214 1.8403];
idx = uint32([1 4 5; 1 3 4; 1 2 3; 1 2 5; 4 5 6; 3 4 6; 2 3 6; 2 5 6]);
areas = [1 1 1 1 1 1 1 1]*pi/2;

[grid,tri] = sphgrid('C1',2);

ok(1) = areequal(phi,grid.phi,1e-8,'abs');
ok(2) = areequal(theta,grid.theta,1e-8,'abs');
ok(3) = areequal(vecs,grid.vecs,1e-8,'abs');
ok(4) = areequal(weights,grid.weights,1e-4,'abs');
ok(5) = all(idx(:)==tri.idx(:));
ok(6) = areequal(areas(:),tri.areas(:),1e-8,'abs');
