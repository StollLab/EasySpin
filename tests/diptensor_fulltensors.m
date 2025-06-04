function ok = test()

rng(2901);

g1 = rand(3,3);
g2 = rand(3,3);
r = [4; 7; -3]/10;

T_test = diptensor(g1,g2,r);

n = r/norm(r);

d12 = 3*(n*n.') - eye(3);
T_ref = -mu0/(4*pi)*bmagn^2/norm(r*1e-9)^3*g1.'*d12*g2/planck/1e6;

ok = areequal(T_test,T_ref,1e-13,'abs');
