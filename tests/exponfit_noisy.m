function ok = test()

% Simple fits to noisy data
%======================================================
rng(45454);
x = linspace(0,10,401);
for q = 1:5
  k = rand*2;
  c = rand*10;
  if rand>0.7; c = -c; end
  y = c*exp(-k*x);
  y = y + 0.05*(max(y)-min(y))*rand(size(x));
  [kf,cf,yf] = exponfit(x,y);
  ok = areequal(kf,k,1e-1,'abs');
end
