function ok = test() %#ok<FNDEF>

% Simulation frame equal to the Zeeman frequency: the Zeeman frequency in
% the simulation frame is zero, so an Sx coherence must not precess
% (single orientation)

Sys.S = 1/2;
Sys.initState = sop(1/2,'x');

Pulse.Type = 'rectangular';
Pulse.tp = 0.001;   % µs
Pulse.Flip = 0.001; % rad, negligible
Exp.Sequence = {Pulse};
Exp.mwFreq = 9.50;
Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree at 9.5 GHz
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = [0 0.02];  % µs
Exp.DetOperator = {'+1'};
Opt.SimFreq = 9.50;
Opt.IntTimeStep = 1e-4;

[~,y] = saffron(Sys,Exp,Opt);

ok = areequal(y,y(1)*ones(size(y)),1e-6,'abs');
