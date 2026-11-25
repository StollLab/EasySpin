function ok = test() %#ok<FNDEF>

% Sys.ZeemanFreq is no longer supported, and the error names it

Pulse.Type = 'rectangular';
Pulse.tp = 0.016;  % µs
Pulse.Flip = pi/2;  % rad

Exp.Sequence = {Pulse 0.1};
Exp.mwFreq = 9.5;  % GHz
Exp.Field = 340;  % mT
Exp.MolFrame = [0 0 0];  % single orientation

Sys = struct('S',1/2,'ZeemanFreq',9.5);  % GHz
try
  saffron(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'ZeemanFreq');
end
