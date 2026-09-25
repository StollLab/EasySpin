function ok = test()

% Negative magnetic field is rejected

Sys.Nucs = '1H';
Sys.A = 5;
Exp.Field = -340;
Exp.Range = [10 20];

try
  salt(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'Exp.Field cannot be negative');
end
