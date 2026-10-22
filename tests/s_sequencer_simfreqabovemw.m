function ok = test()

% Pulses below the simulation frame (negative frequencies in that frame) need
% the same integration time step as pulses equally far above it.

Pulse.Type = 'rectangular';
Pulse.tp = 0.02;  % µs
Pulse.Flip = pi/2;

Exp.Sequence = {Pulse 0.1};
Exp.Field = 340;
Exp.mwFreq = 9.5;  % GHz

Opt.SimFreq = 10;  % GHz, 0.5 GHz above the pulses
[~,~,OptAbove] = runprivate('s_sequencer',Exp,Opt);
Opt.SimFreq = 9;   % GHz, 0.5 GHz below the pulses
[~,~,OptBelow] = runprivate('s_sequencer',Exp,Opt);

ok(1) = OptAbove.IntTimeStep>0;
ok(2) = OptAbove.IntTimeStep==OptBelow.IntTimeStep;

% A sufficiently small user-provided time step is accepted
Opt.SimFreq = 10;
Opt.IntTimeStep = 1e-5;  % µs
try
  runprivate('s_sequencer',Exp,Opt);
  ok(3) = true;
catch
  ok(3) = false;
end
