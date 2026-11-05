function ok = test() %#ok<FNDEF>

% Without Exp.Field, saffron gives an informative error

Pulse.tp = 0.02;  % µs
Pulse.Flip = pi/2;
Exp.Sequence = {Pulse 0.1};
Exp.mwFreq = 9.5;  % GHz

Sys.S = 1/2;
try
  saffron(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'Exp.Field is missing');
end
