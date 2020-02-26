function ok = test()

% Simple fits to noisefree data

x = linspace(0,10,201);
rng(4);
for q = 1:5
  k = rand*2;
  c = rand*10;
  if rand>0.7; c = -c; end
  y = c*exp(-k*x);
  kf = exponfit(x,y);
  ok(q) = areequal(kf,k,1e-6,'abs');
end
