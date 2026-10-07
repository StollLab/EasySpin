function ok = test()

% Sys.ZeemanFreq overrides Sys.g for all g input forms

Pulse.Type = 'rectangular';
Pulse.tp = 0.01; % µs
Pulse.Flip = pi/2;

Exp.Sequence = {Pulse 0.05};
Exp.Field = 1240;
Exp.TimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;
Exp.DetSequence = [0 1];

Sys.S = 1/2;
Sys.ZeemanFreq = 33.5;
[~,signal0] = spidyan(Sys,Exp);

gList = {[2.0 2.1], [2.0 2.1 2.2], diag([2.0 2.1 2.2]), [2.0 2.1 2.2 0.01 0.02 0.03]};
for k = 1:numel(gList)
  Sys.g = gList{k};
  [~,signal] = spidyan(Sys,Exp);
  ok(k) = areequal(signal,signal0,1e-10,'abs');
end
