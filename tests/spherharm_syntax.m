function ok = test()

q = 1000;
theta = rand(1,q);
phi = rand(1,q);
L = 5;
M = 4;

y = spherharm(L,M,theta,phi);
y = spherharm(L,M,theta,phi,'r');

ok = true;
