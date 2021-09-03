function ok = test()

N = 100;
nOctants = 8;  % full sphere
v = sphrand(N,nOctants);

ok(1) = size(v,1)==3 && size(v,2)==N;
ok(2) = min(v(3,:))<0;  % check for negative z coordinates
ok(3) = min(v(1,:))<0;  % check for negative x coordinates
ok(4) = min(v(2,:))<0;  % check for negative y coordinates
