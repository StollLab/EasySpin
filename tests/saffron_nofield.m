function ok = test() %#ok<FNDEF>

% Without Exp.Field, saffron gives an informative error

p90.tp = 0.02;  % µs
p90.Flip = pi/2;  % rad
tau = 0.1;  % µs
Exp.Sequence = {p90 tau};
Exp.mwFreq = 9.5;  % GHz

Sys.S = 1/2;
try
  saffron(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'Exp.Field is missing');
end
