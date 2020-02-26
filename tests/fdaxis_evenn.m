function ok = test()

% Test for even N
%------------------------------------------------
dT = 10;

N = 532;
fx1 = fdaxis(dT,N);
fx2 = linspace(-N/2,N/2-1,N)*(1/N/dT);

ok = areequal(fx1,fx2,1e-10,'abs');
