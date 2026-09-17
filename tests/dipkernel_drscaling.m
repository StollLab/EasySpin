function ok = test()

% The kernel is scaled by the distance-grid spacing dr. Halving dr
% (while keeping the same shared grid points) must halve K at those
% shared points.

t = 0.5;  % us

r1 = 2:0.1:4;    % dr = 0.1 nm
r2 = 2:0.05:4;   % dr = 0.05 nm, r2(1:2:end) coincides with r1

K1 = dipkernel(t,r1);
K2 = dipkernel(t,r2);
K2_shared = K2(1:2:end);

ok = areequal(K2_shared,0.5*K1,1e-10,'rel');

end
