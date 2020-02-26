function ok = test()

% Hamming window

N = 2345;
w0 = apowin('ham',N);
x = linspace(-1,1,N);
w1 = 0.54 + 0.46*cos(pi*x(:));

ok = areequal(w0,w1,1e-10,'abs');
