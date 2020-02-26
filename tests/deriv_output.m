function ok = test()

N = 2345;
y = rand(1,N);
x = 1:N;
a = deriv(x,y);

ok = numel(a)==numel(y);

