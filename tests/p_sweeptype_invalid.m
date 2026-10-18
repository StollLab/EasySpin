function ok = test()

% Invalid or inconsistent mw frequency, field and sweep range fields are rejected

n = 1024;
Exps = {
  struct('nPoints',n), 'or Exp.Field (for a frequency sweep).'  % neither given
  struct('nPoints',n,'mwFreq',9.5,'Field',340), 'but not both'
  struct('nPoints',n,'mwFreq',9.5,'mwRange',[9 10]), 'Exp.mwRange cannot be used in a field sweep'
  struct('nPoints',n,'mwFreq',9.5,'mwCenterSweep',[9.5 1]), 'Exp.mwCenterSweep cannot be used in a field sweep'
  struct('nPoints',n,'Field',340,'Range',[300 400]), 'Exp.Range cannot be used in a frequency sweep'
  struct('nPoints',n,'Field',340,'CenterSweep',[340 10]), 'Exp.CenterSweep cannot be used in a frequency sweep'
  struct('nPoints',n,'mwFreq',[9 10]), 'Uninterpretable microwave frequency'
  struct('nPoints',n,'mwFreq',-9.5), 'Uninterpretable microwave frequency'
  struct('nPoints',n,'Field',NaN), 'Uninterpretable magnetic field'
  struct('nPoints',n,'Field',-340), 'Exp.Field cannot be negative'
  struct('nPoints',n,'Field',0), 'Exp.Field cannot be zero'
  struct('nPoints',1,'mwFreq',9.5), 'Problem with Exp.nPoints'
  struct('nPoints',100.5,'mwFreq',9.5), 'Problem with Exp.nPoints'
  struct('nPoints',NaN,'mwFreq',9.5), 'Problem with Exp.nPoints'
  struct('nPoints',Inf,'mwFreq',9.5), 'Problem with Exp.nPoints'
  struct('nPoints',n,'mwFreq',9.5,'Mode','parallel'), 'Exp.Mode is no longer supported'
  };

allowZeroField = false;
for k = size(Exps,1):-1:1
  try
    runprivate('p_sweeptype',Exps{k,1},allowZeroField);
    ok(k) = false;
  catch err
    ok(k) = contains(err.message,Exps{k,2});
  end
end

end
