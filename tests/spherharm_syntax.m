function ok = test()

q = 1000;
theta = rand(1,q);
phi = rand(1,q);
y = spherharm(5,4,theta,phi);

ok = true;
