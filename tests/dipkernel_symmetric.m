function ok = test()

% dipkernel uses abs(t) internally, so the kernel must be an even
% function of t.

t = linspace(-2,2,21);  % us
r = 3.5;  % nm

K1 = dipkernel(t,r);
K2 = dipkernel(-t,r);

ok = areequal(K1,K2,1e-12,'abs');

end
