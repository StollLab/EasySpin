function ok = test()

% Zero relaxation times are ignored

Sys.S = 1/2;

Pulse.Type = 'rectangular';
Pulse.tp = 0.02;
Pulse.Flip = pi/2;

Exp.Sequence = {Pulse 0.1};
Exp.mwFreq = 9.5;
Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree
Exp.DetSequence = [0 1];

Opt.SimFreq = 9.4;

[~,y0] = spidyan(Sys,Exp,Opt);

Sys.T1 = 0;
Sys.T2 = 0;
[~,y] = spidyan(Sys,Exp,Opt);

ok = areequal(y,y0,1e-10,'rel');
