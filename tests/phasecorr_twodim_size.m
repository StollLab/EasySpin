function ok = test()

% For 2D input, the output size must match the input size exactly, and
% the phase found must minimize the global (not per-row) imaginary
% energy, matching a brute-force grid search.

t = linspace(0,3,100);  % µs
f = 2.445;  % MHz
V = cos(2*pi*f*t).*exp(1i*deg2rad(12));
rng(5);
V = V + 0.1*randn(5,numel(t));

[Vph,ph] = phasecorr(V);
ok(1) = isequal(size(Vph),size(V));
ok(2) = areequal(Vph,V.*exp(1i*ph),1e-12,'abs');

grid = linspace(-pi/2,pi/2,401);
cost = zeros(size(grid));
Vflat = V(:);
for k = 1:numel(grid)
  cost(k) = sum(imag(Vflat*exp(1i*grid(k))).^2);
end
[~,idx] = min(cost);
ok(3) = areequal(ph,grid(idx),2e-2,'abs');

end
