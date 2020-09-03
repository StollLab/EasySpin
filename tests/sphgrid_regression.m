function [ok,data] = test(opt,olddata)

% Regression test

Groups = {'D2h','C2h','Dinfh','Ci'};
nKnots = [30 10 40 12];
for g = 1:numel(Groups)
  grid = sphgrid(Groups{g},nKnots(g));
  v{g} = grid.vecs;
end

data.v = v;

if isempty(olddata)

  ok = [];

else

  ok = true;
  for k=1:numel(v)
    if any(abs(v{k}(:)-olddata.v{k}(:))>1e-10)
      ok = false;
      break
    end
  end

end
