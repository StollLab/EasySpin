function ok = test()

% At t=0, the dipolar kernel value is 1 (before the dr distance-grid
% scaling factor is applied).

r_single = 3.5;  % nm, single distance -> no dr scaling applied
K0 = dipkernel(0,r_single);
ok(1) = areequal(K0,1,1e-12,'abs');

r = 2:0.1:4;  % nm, equidistant grid -> dr scaling applied
dr = r(2)-r(1);
K0r = dipkernel(0,r);
ok(2) = areequal(K0r,dr*ones(size(r)),1e-12,'abs');

end
