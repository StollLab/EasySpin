function [err,data] = test(opt,olddata)

% Test whether associated Legendre polynomials are orthogonal
%============================================================

thr = 1e-8;
v = quad(@integrand1,-1,1,thr);
n = 6; m = 2;
v0 = 2/(2*n+1) * factorial(n+m)/factorial(n-m);
err = abs(v-v0)>10*thr;

data = [];

function z = integrand1(z)
z = plegendre(6,2,z).*plegendre(6,2,z);
