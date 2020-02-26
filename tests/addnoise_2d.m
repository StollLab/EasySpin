function ok = test()

%======================================================
% Check whether returned 2D array size is correct
%======================================================

rng(1);
y = rand(14,1653);

SNR = 20;
noisemodel = {'f','u','n'};

for k = 1:3
  y_ = addnoise(y,SNR,noisemodel{k});
  ok(k) = all(size(y_)==size(y));
end
