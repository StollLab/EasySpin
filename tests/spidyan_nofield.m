function ok = test()

% Without Exp.Field, spidyan needs Sys.ZeemanFreq. With neither Sys.g nor
% Sys.ZeemanFreq, it must issue an informative error.

Pulse.tp = 0.02;  % µs
Pulse.Flip = pi/2;
Exp.Sequence = {Pulse 0.1};
Exp.mwFreq = 9.5;  % GHz

Sys.S = 1/2;
try
  spidyan(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'Exp.Field is required');
end
