function [err,data] = test(opt,olddata)

% Test 2: simple fits to noisy data
%======================================================
x = linspace(0,10,401);
ok = 1;
for q = 1:50
  k = rand*2;
  c = rand*10;
  if rand>0.7; c = -c; end
  y = c*exp(-k*x);
  y = y + 0.05*(max(y)-min(y))*rand(size(x));
  [kf,cf,yf] = exponfit(x,y);
  if abs(kf-k)>1e-1
    ok = 0;
    break;
  end
end
err = ~ok;
data = [];
