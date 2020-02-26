function ok = test()

% Test for odd N
%------------------------------------------------
dT = rand;

N = 345;
fx1 = fdaxis(dT,N);
fx2 = (-floor(N/2):floor(N/2))*(1/N/dT);

ok = areequal(fx1,fx2,1e-10,'abs');
