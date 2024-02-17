function ok = test()

% One 3x3 g tensor and and one isotropic g value

rng(723)

gtensor = diag([2 2.1 2.2]) + rand(3,3)*0.1;
giso = gfree;

r = [-1;6;2]/10;  % nm

T_test = diptensor(gtensor,giso,r);


n = r/norm(r);
d12 = 3*(n*n.') - eye(3);
T_ref = -mu0/(4*pi)*bmagn^2/norm(r*1e-9)^3*gtensor.'*d12*giso/planck/1e6;

ok = areequal(T_test,T_ref,1e-12,'abs');
