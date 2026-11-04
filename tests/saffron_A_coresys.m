function ok = test() %#ok<FNDEF>

% Check for a fixed bug, where saffron crashed if Sys.A_ was given and the
% excitation window of ideal pulses (Opt.IdealPulses) triggered the removal
% of a nucleus for coreSys.

Sys.S = 1/2;
Sys.Nucs = '1H,1H';
Sys.lwEndor = 0.1;

Exp.Field = 350;
Exp.Sequence = 'MimsENDOR';

Sys.A_ = [200 20; 3 2];

Exp.tau = 444;
Exp.Range = [13 16.4];
Exp.mwFreq = 9.7;
Exp.tp = 0.02;  % µs, excitation window 1/tp = 50 MHz
Opt.IdealPulses = true;

[x,y] = saffron(Sys,Exp,Opt);

ok = true;
