function ok = test()

N = 1000;
x = 1:N;
dx = x(2) - x(1);
y = rand(1,N);
dydx = diff(y)/dx;
dydx = (dydx([1 1:end]) + dydx([1:end end]))/2;
a = deriv(x,y);
ok = areequal(dydx,a,1e-10,'abs');

