function ok = test()

% Test whether associated Legendre polynomials are orthogonal

thr = 1e-8;

integrand1 = @(z) plegendre(6,2,z).*plegendre(6,2,z);
v = integral(integrand1,-1,1,'AbsTol',thr);

L = 6; M = 2;
v0 = 2/(2*L+1) * factorial(L+M)/factorial(L-M);

ok = areequal(v,v0,10*thr,'abs');

