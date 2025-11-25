function ok = test()

ph = 0:15:360;  % deg
t = linspace(0,3,1001);  % µs
f = 2.445;  % MHz
V = cos(2*pi*f*t);

for k = 1:numel(ph)
  V_ = V.*exp(1i*deg2rad(ph(k)));
  [~,phfit(k)] = phasecorr(V_);
end

imagerr = imag(exp(1i*deg2rad(ph)).*exp(1i*phfit));

ok = abs(imagerr)<1e-3;

end
