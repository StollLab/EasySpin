function [err,data] = test(opt,olddata)

% Syntax check

Sy = {'D2h','C2h','Ci','Dinfh'};
nK = 10;
for k=1:numel(Sy)
  v = sphgrid(Sy{k},nK);
  [v,w] = sphgrid(Sy{k},nK,'c');
  [p,t] = sphgrid(Sy{k},nK);
  [p,t,w] = sphgrid(Sy{k},nK);
end

err = 0;
data = [];
