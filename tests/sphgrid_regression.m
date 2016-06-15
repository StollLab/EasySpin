function [err,data] = test(opt,olddata)

% Regression test

Groups = {'D2h','C2h','Dinfh','Ci'};
nKnots = [30 10 40 12];
for g = 1:numel(Groups)
  v{g} = sphgrid(Groups{g},nKnots(g));
end

data.v = v;

if isempty(olddata)

  err = [];

else

  err = 0;
  for k=1:numel(v)
    if any(abs(v{k}(:)-olddata.v{k}(:))>1e-10)
      err = 1;
      break;
    end
  end

end
