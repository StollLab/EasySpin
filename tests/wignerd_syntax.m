function ok = test()

J = 15/2;
m1 = -7/2;
m2 = 3/2;
Jm1m2 = [J m1 m2];

rng(1);

alpha = rand*2*pi;
beta = rand*pi;
gamma = rand*2*pi;

wignerd(J,alpha,beta,gamma);
wignerd(J,alpha,beta,gamma,'+');
wignerd(J,alpha,beta,gamma,'-');
wignerd(J,beta);
wignerd(J,beta,'+');
wignerd(J,beta,'-');
wignerd(Jm1m2,alpha,beta,gamma);
wignerd(Jm1m2,alpha,beta,gamma,'+');
wignerd(Jm1m2,alpha,beta,gamma,'-');
wignerd(Jm1m2,beta);
wignerd(Jm1m2,beta,'+');
wignerd(Jm1m2,beta,'-');

ok = true;
