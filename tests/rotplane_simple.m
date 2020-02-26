function ok = test()

chi = linspace(0,2*pi,10);
n = [0;0;1];

v1 = rotplane(n,chi);
v2 = ang2vec(-chi,pi/2);

ok = areequal(v1,v2,1e-10,'abs');
