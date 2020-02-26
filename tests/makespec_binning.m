function ok = test()

% Correct binning

ok = true;
Amp = 4;
for k = 1:100
  Pos = rand;
  [x,s] = makespec([0 1],101,Pos,Amp);
  idx = find(s>0);
  dist = abs(s-Pos);
  if abs(dist(idx)-min(dist))<1e-10
    ok = false;
    break
  end
  if s(idx)~=Amp
    ok = false;
    break
  end
end

