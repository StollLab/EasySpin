function ok = test()

% Negative magnetic field is rejected

Sys.g = 2;
Exp.Field = -340;
Exp.Range = [9 10];

try
  resfreqs_perturb(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'Exp.Field cannot be negative');
end
