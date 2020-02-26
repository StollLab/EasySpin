function ok = test()

% alpha gamma degeneracy

phi = pi/7;
R0 = erot(phi,0,0);
R1 = erot(0,0,phi);
R2 = erot(phi/2,0,phi/2);

thr = 1e-10;
ok = areequal(R0,R1,thr,'abs') && areequal(R0,R2,thr,'abs');
