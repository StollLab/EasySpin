function ok = test()

N = 654;
y = rand(1,N);
x = 1:N;
a = deriv(y);
b = deriv(x,y);

ok = areequal(a,b,1e-10,'abs');
