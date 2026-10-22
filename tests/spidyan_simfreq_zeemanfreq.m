function ok = test()

% Simulation frame equal to Sys.ZeemanFreq: the Zeeman frequency in the
% simulation frame is zero, so an Sx coherence must not precess.
% (Previously, the frame shift was subtracted from Sys.ZeemanFreq before
% its conversion to g, which skips zero entries, so the spin kept the
% wrong g value.)

Sys.S = 1/2;
Sys.ZeemanFreq = 9.50;  % GHz
Sys.initState = sop(1/2,'x');

Pulse.Type = 'rectangular';
Pulse.tp = 0.001;   % µs
Pulse.Flip = 0.001; % rad, negligible
Exp.Sequence = {Pulse 0.02};
Exp.Field = 340;
Exp.mwFreq = 9.50;
Exp.DetSequence = [0 1];
Exp.DetOperator = {'+1'};
Opt.SimFreq = 9.50;
Opt.IntTimeStep = 1e-4;

[~,y] = spidyan(Sys,Exp,Opt);

ok = areequal(y,y(1)*ones(size(y)),1e-6,'abs');
