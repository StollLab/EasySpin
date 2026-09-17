function ok = test()

% The dipolar frequency, and therefore the phase entering the kernel,
% scales linearly with the product gA*gB. Scaling gA by a factor f
% (with gB fixed) must be equivalent to scaling t by the same factor f
% with the default g values.

t = linspace(0.05,2,15);  % us
r = 2:0.2:4;  % nm
f = 1.7;

K1 = dipkernel(t,r,[f*gfree gfree]);
K2 = dipkernel(f*t,r);

ok = areequal(K1,K2,1e-10,'rel');

end
