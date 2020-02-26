function ok = test()

% Syntax check

SymmGroup = {'D2h','C2h','Ci','Dinfh'};
nKnots = 10;
for g = 1:numel(SymmGroup)
  v = sphgrid(SymmGroup{g},nKnots);
  [v,w] = sphgrid(SymmGroup{g},nKnots,'c');
  [p,t] = sphgrid(SymmGroup{g},nKnots);
  [p,t,w] = sphgrid(SymmGroup{g},nKnots);
end

ok = true;
