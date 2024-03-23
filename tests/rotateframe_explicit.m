function ok = test()

angles0 = pi*[0.4 0.76 -0.98];  % Euler angles, radians
nRot = [1 2 3];  % rotation axis
rho = -3*pi/7;  % rotation angle, radians

% Explicit rotation
R = rotaxi2mat(nRot,rho);
xyz = erot(angles0).';
xyz_rotated = R.'*xyz;
nocheck = true;
angles_ref = eulang(xyz_rotated.',nocheck);

% Rotation using rotateframe()
angles = rotateframe(angles0,nRot,rho);

ok = areequal(angles_ref,angles,1e-8,'abs');
