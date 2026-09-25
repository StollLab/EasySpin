function ok = test()

% Negative magnetic field is rejected

Sys.Nucs = '1H';
Sys.A = 5;
Exp.Field = -340;
Exp.Sequence = '2pESEEM';
Exp.dt = 0.01;

try
  saffron(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'Exp.Field cannot be negative');
end
