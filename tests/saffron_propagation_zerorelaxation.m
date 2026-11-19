function ok = test() %#ok<FNDEF>

% Zero relaxation times are ignored (single orientation, FID after a pi/2
% pulse)

Sys.S = 1/2;

Pulse.Type = 'rectangular';
Pulse.tp = 0.02;
Pulse.Flip = pi/2;

Exp.Sequence = {Pulse};
Exp.mwFreq = 9.5;
Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = [0 0.1];  % µs

Opt.SimFreq = 9.4;

[~,y0] = saffron(Sys,Exp,Opt);

Sys.T1 = 0;
Sys.T2 = 0;
[~,y] = saffron(Sys,Exp,Opt);

ok = areequal(y,y0,1e-10,'rel');
