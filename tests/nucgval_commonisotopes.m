function ok = test(opt)

% Test gn values of common isotopes

Isotopes = {'1H','2H','13C','14N','15N','33S','59Co','63Cu'};
gn = [5.58569468,0.8574382,1.4048236,0.40376100,-0.56637768,0.429214,1.322,1.4824];

for k = 1:numel(Isotopes)
  ngn(k) = nucgval(Isotopes{k});
end
if opt.Display
  ngn-gn
end

ok = areequal(ngn,gn,1e-4,'rel');
