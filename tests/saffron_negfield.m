function ok = test()

% Negative magnetic field is rejected

Sys.Nucs = '1H';
Sys.A = 5;  % MHz
Exp.Field = -340;  % mT
Exp.Sequence = '2pESEEM';
Exp.dt = 0.01;  % µs

try
  saffron(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'Exp.Field cannot be negative');
end
