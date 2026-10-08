function ok = test()

% Opt.Freq2Field is rejected

Sys.g = 2;
Sys.lw = 1;
Exp.mwFreq = 9.5;
Exp.Range = [330 350];
Opt.Freq2Field = 0;

try
  pepper(Sys,Exp,Opt);
  ok = false;
catch err
  ok = contains(err.message,'Freq2Field');
end
