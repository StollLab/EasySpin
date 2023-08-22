function ok = test()

Frame = pi*[0.4 0.76 -0.98];  % Euler angles, radians
nRot = [1 2 3];  % rotation axis
rho = -3*pi/7;  % rotation angle, radians

% Explicit rotation
R = rotaxi2mat(nRot,rho);
xyz = erot(Frame).';
xyz_rotated = R.'*xyz;
skipFitting = true;
angles0 = eulang(xyz_rotated.',skipFitting);

% Rotation using rotateframe()
angles1 = rotateframe(Frame,nRot,rho);

ok = areequal(angles0,angles1,1e-8,'abs');
