function ok = test()

% Test against explicit values

% Cartesian rank-2 tensor
T = [1 2 3; 4 5 6; 7 8 9];

% Calculate spherical tensors
[T0,T1,T2] = tensor_cart2sph(T);

% Reference spherical tensors
T0ref = -sqrt(75);
T1ref = [-2-1i; sqrt(2)*1i; -2+1i];
T2ref = [-2+3i; -5-7i; sqrt(24); 5-7i; -2-3i];

thr = 1e-14;
ok(1) = areequal(T0,T0ref,thr,'rel');
ok(2) = areequal(T1,T1ref,thr,'rel');
ok(3) = areequal(T2,T2ref,thr,'rel');
