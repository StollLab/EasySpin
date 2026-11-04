function ok = test()

% Simulation frame equal to the Zeeman frequency: the Zeeman frequency in
% the simulation frame is zero, so an Sx coherence must not precess.

Sys.S = 1/2;
Sys.initState = sop(1/2,'x');

Pulse.Type = 'rectangular';
Pulse.tp = 0.001;   % µs
Pulse.Flip = 0.001; % rad, negligible
Exp.Sequence = {Pulse 0.02};
Exp.mwFreq = 9.50;
Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree at 9.5 GHz
Exp.DetSequence = [0 1];
Exp.DetOperator = {'+1'};
Opt.SimFreq = 9.50;
Opt.IntTimeStep = 1e-4;

[~,y] = spidyan(Sys,Exp,Opt);

ok = areequal(y,y(1)*ones(size(y)),1e-6,'abs');
