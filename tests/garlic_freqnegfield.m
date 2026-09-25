function ok = test()

% Negative magnetic field is rejected

Sys.g = 2;
Sys.lw = 1;
Exp.Field = -340;
Exp.mwRange = [9 10];

try
  garlic(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'Exp.Field cannot be negative');
end
