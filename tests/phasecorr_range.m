function ok = test()

% The returned phase must always lie within the documented range
% (-pi/2,pi/2], for both ignoreOffset settings, across a variety of
% random signals.

rng(1);
nTrials = 30;
ok = true(1,nTrials);
for k = 1:nTrials
  n = randi([10 200]);
  V = randn(1,n) + 1i*randn(1,n) + (rand-0.5)*10i;
  [~,ph1] = phasecorr(V,false);
  [~,ph2] = phasecorr(V,true);
  inrange1 = ph1>-pi/2-1e-10 && ph1<=pi/2+1e-10;
  inrange2 = ph2>-pi/2-1e-10 && ph2<=pi/2+1e-10;
  ok(k) = inrange1 && inrange2;
end

end
