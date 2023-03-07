function ok = test()

% 120 degree rotation of z around around 111 should give x

n = [1;1;1];
rho = 2*pi/3;
R = rotaxi2mat(n,rho);

z = [0;0;1];
x = [1;0;0];
z_rot = R.'*z;  % active rotation

ok = areequal(z_rot,x,1e-7,'abs');
