function ok = test()

Frame = pi*[0.4 0.76 -0.98];  % Euler angles, radians
nRot = [1 2 3];  % rotation axis
rho = -3*pi/7;  % rotation angle, radians

% Explicit rotation
R = rotaxi2mat(nRot,rho);
xyzC_L = erot(Frame);
xyzC_L_rotated = R*xyzC_L;
skipFitting = true;
angles0 = eulang(xyzC_L_rotated,skipFitting);

% Rotation using rotateframe()
angles1 = rotateframe(Frame,nRot,rho);

ok = areequal(angles0,angles1,1e-8,'abs');
