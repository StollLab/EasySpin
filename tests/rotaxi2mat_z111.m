function ok = test()

% 120 degree rotation of z around around 111 should give y

n = [1;1;1];
rho = 2*pi/3;
R = rotaxi2mat(n,rho);

z = [0;0;1];
y = [0;1;0];

ok = areequal(R*z,y,1e-7,'abs');
