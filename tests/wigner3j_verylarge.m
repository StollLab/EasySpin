function ok = test()

% Test for very large parameters

a(1) = wigner3j(5000,4000,1000,25,-15,-10);
b(1) = -0.0013903154820072137512;
a(2) = wigner3j(500,1200,1500,200,-150,-50);
b(2) = 0.00060418496059993375494;

ok = areequal(a,b,1e-9,'abs');
