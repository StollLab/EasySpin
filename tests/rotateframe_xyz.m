function ok = test()

% Test small rotations of an aligned frame around x, y anx z axes.

angles0 = [0 0 0];  % Euler angles, radians
rho = 0.2;         % rotation angle, radians

angles = rotateframe(angles0,'y',rho);
angles_ref = [0 rho 0];
ok(1) = areequal(angles,angles_ref,1e-8,'abs');

angles = rotateframe(angles0,'z',rho);
angles_ref = [rho 0 0];
ok(2) = areequal(angles,angles_ref,1e-8,'abs');

angles = rotateframe(angles0,'x',rho);
angles_ref = [2*pi-pi/2 rho pi/2];
ok(3) = areequal(angles,angles_ref,1e-8,'abs');
