function ok = test()

% Interface test
dT = 0.45;
N = 200;
tx = (0:N-1)*dT;
fx1 = fdaxis(dT,N);
fx2 = fdaxis(tx);

ok = areequal(fx1,fx2,1e-10,'abs');
