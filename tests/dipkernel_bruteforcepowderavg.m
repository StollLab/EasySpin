function ok = test()

% Cross-validate the closed-form Fresnel-integral expression against a
% direct numerical powder average of the dipolar oscillation
% cos(wdd*t*(1-3*cos(theta)^2)) over cos(theta) in [0,1].

r = 3.5;  % nm, single distance -> no dr scaling
tvec = [0.1 0.3 0.7 1.5];  % µs

gA = gfree;
gB = gfree;
D = (mu0/4/pi)*bmagn^2*gA*gB/planck*1e21;  % MHz nm^3
wdd = 2*pi*D/r^3;  % rad/µs

K_ref = zeros(size(tvec));
for k = 1:numel(tvec)
  t = tvec(k);
  K_ref(k) = integral(@(x) cos(wdd*t*(1-3*x.^2)),0,1);
end

K = dipkernel(tvec,r);

ok = areequal(K(:),K_ref(:),1e-6,'abs');

end
