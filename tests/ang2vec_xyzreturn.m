function ok = test()

% Expansion of a single angle to a vector

N = 100;
phi = rand(1,N);
theta = rand(1,N);

V = ang2vec(phi,theta);
[x,y,z] = ang2vec(phi,theta);
V1 = [x;y;z];

ok = areequal(V,V1);
