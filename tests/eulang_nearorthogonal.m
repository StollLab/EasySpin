function ok = eulang_nearorthogonal()

% Test whether eulang works with a near-orthogonal matrix

rng(4545);

% Generate orthogonal matrix
angles0 = [10 40 80]*pi/180;
R = erot(angles0);

% Add noise
R = R + 0.001*rand(3,3);

% Calculate Euler angles
angles = eulang(R);

ok = areequal(angles,angles0,5e-3,'abs');
