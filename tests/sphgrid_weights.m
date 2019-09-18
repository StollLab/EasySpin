function [err,data] = test(opt,olddata)

% Verify that weights always add up to 4*pi

nKnots = 10;

SymmGroups = {'C1','Ci','C2h','S6','C4h','C6h','D2h','Th',...
  'D3d','D4h','Oh','D6h','Dinfh','O3'};

for g = 1:numel(SymmGroups)
  [Vecs,w_] = sphgrid(SymmGroups{g},nKnots,'c');
  wsum(g) = sum(w_);
end
wsum = wsum/(4*pi);

err = any(abs(wsum-1)>1e-10);

data = [];
