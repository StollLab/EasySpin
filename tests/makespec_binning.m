function [err,data] = test(opt,olddata)

% Test 2: correct binning
%====================================================
err = 0;
Amp = 4;
for k=1:100
  Pos = rand;
  [x,s] = makespec([0 1],101,Pos,Amp);
  idx = find(s>0);
  dist = abs(s-Pos);
  if abs(dist(idx)-min(dist))<1e-10
    err(2) = 1;
    break;
  end
  if s(idx)~=Amp
    err = 1;
    break;
  end
end
data = [];
