function ok = test()

% Negative sweep ranges are accepted or rejected depending on allowNegative

Exp.CenterSweep = [0 100];
R = runprivate('p_sweeprange',Exp,false,true);
ok = isequal(R,[-50 50]);

try
  runprivate('p_sweeprange',Exp,false,false);
  ok = false;
catch err
  ok = ok && contains(err.message,'Sweep range cannot be negative! Check Exp.CenterSweep or Exp.Range.');
end
