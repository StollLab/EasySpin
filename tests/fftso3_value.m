function ok = test()

% Define function
f = @(a,b,c)exp(-0.2*wignerd([2 1 1],a,b,c)-0.4*wignerd([1 0 -1],a,b,c));

% Calculate FFT
Lmax = 12;
c = fftso3(f,Lmax);

% Pick random orientation
rng(252);
alpha = rand*2*pi;
beta = rand*pi;
gamma = rand*2*pi;

% Evaluate function
y0 = f(alpha,beta,gamma);

% Evaluate function using FFT coefficients
y = 0;
for L = 0:Lmax
  D = wignerd(L,alpha,beta,gamma);
  y = y + sum(sum(c{L+1}.*D(end:-1:1,end:-1:1)));
end

% Comparison
ok = areequal(y,y0,1e-8,'abs');
