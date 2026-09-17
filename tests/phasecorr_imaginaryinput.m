function ok = test()

% A purely imaginary signal (no offset) must be rotated fully onto the
% real axis.

t = linspace(0,3,1001);  % µs
f = 2.445;  % MHz
x = cos(2*pi*f*t);
V = 1i*x;

[Vph,ph] = phasecorr(V);
ok(1) = areequal(abs(ph),pi/2,1e-8,'abs');
ok(2) = areequal(imag(Vph),zeros(size(V)),1e-8,'abs');
ok(3) = areequal(abs(real(Vph)),abs(x),1e-8,'abs');

end
