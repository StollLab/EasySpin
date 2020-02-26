function ok = test()

% Weights must add up to 4*pi

Sy = {'D2h','C2h','Ci','Dinfh'};
nK = 10;
for k=1:numel(Sy)
  [dum,dum,w] = sphgrid(Sy{k},nK);
  Integral(k) = sum(w);
end

ok = all(abs(Integral/(4*pi)-1)<1e-4);
