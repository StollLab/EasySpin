function ok = test()

% Explicitly check very small grid C2h

phi = [0 0 0.5 0 0.25 0.5 0.75]*pi;
theta = [0 0.25 0.25 0.5 0.5 0.5 0.5]*pi;
c = 1/sqrt(2);
vecs = [0 0 1; c 0 c; 0 c c; 1 0 0; c c 0; 0 1 0; -c c 0].';
weights = [0.9566, 3.4004, 3.4004, 1.2022, 1.2022, 1.2022, 1.2022];
triidx = uint32([1 2 3; 1 2 3; 2 3 5; 2 3 7; 2 4 5; 2 4 7; 3 5 6; 3 6 7]);
triareas = [1.3593 1.3593 2.2051 2.2051 1.3593 1.3593 1.3593 1.3593];

[grid,tri] = sphgrid('C2h',3);

ok(1) = areequal(phi,grid.phi,1e-8,'abs');
ok(2) = areequal(theta,grid.theta,1e-8,'abs');
ok(3) = areequal(vecs,grid.vecs,1e-8,'abs');
ok(4) = areequal(weights,grid.weights,1e-4,'abs');
ok(5) = all(triidx(:)==tri.idx(:));
ok(6) = areequal(triareas(:),tri.areas(:),1e-4,'abs');
