function ok = test()

% Cross-validate the analytic phase against a brute-force grid search
% that minimizes sum(imag(V*exp(1i*phi)).^2) independently, for both
% ignoreOffset settings.

t = linspace(0,3,100);  % µs
f = 2.445;  % MHz
V = cos(2*pi*f*t).*exp(1i*deg2rad(23)) + 0.3i;

grid = linspace(-pi/2,pi/2,401);
tol = 2e-2;  % a few times the grid spacing

[~,ph1] = phasecorr(V,false);
cost1 = zeros(size(grid));
for k = 1:numel(grid)
  cost1(k) = sum(imag(V*exp(1i*grid(k))).^2);
end
[~,idx1] = min(cost1);
ok(1) = areequal(ph1,grid(idx1),tol,'abs');

Vzm = V - mean(V);
[~,ph2] = phasecorr(V,true);
cost2 = zeros(size(grid));
for k = 1:numel(grid)
  cost2(k) = sum(imag(Vzm*exp(1i*grid(k))).^2);
end
[~,idx2] = min(cost2);
ok(2) = areequal(ph2,grid(idx2),tol,'abs');

end
