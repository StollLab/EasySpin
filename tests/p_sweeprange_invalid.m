function ok = test()

% Invalid sweep ranges are rejected

Exps = {struct('CenterSweep',[340 -10]), ...
        struct('CenterSweep',[340 0]), ...
        struct('CenterSweep',340), ...
        struct('CenterSweep',[340 10 1]), ...
        struct('CenterSweep',[NaN 10]), ...
        struct('Range',[400 300]), ...
        struct('Range',[300 Inf])};

ok = true;
for k = 1:numel(Exps)
  try
    runprivate('p_sweeprange',Exps{k},false,true);
    ok = false;
  catch err
    ok = ok && contains(err.message,'Invalid sweep range! Check Exp.CenterSweep or Exp.Range.');
  end
end
