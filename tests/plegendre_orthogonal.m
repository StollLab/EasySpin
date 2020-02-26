function ok = test()

% Test whether associated Legendre polynomials are orthogonal

thr = 1e-8;

integrand1 = @(z) plegendre(6,2,z).*plegendre(6,2,z);
v = integral(integrand1,-1,1,'AbsTol',thr);
n = 6; m = 2;
v0 = 2/(2*n+1) * factorial(n+m)/factorial(n-m);

ok = abs(v-v0)<10*thr;

