function ok = test()

% System ------------------------
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.tp = 0.03;
Pulse.Flip = pi;

Exp.Sequence = {Pulse};
Exp.mwFreq = 33.5;
Exp.DetSequence = [1]; 

Exp.DetOperator = 'z1';
[sig1] = spidyan(Sys,Exp);

Exp.DetOperator = {'z1'};
[sig2] = spidyan(Sys,Exp);

Exp.DetOperator = [1/2 0; 0 -1/2];
[sig3] = spidyan(Sys,Exp);

ok = isequal(sig1,sig2) && isequal(sig1,sig3);
