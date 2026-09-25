function ok = test()

% Resonance fields at negative fields are the mirror images of those at
% positive fields

Sys.g = [2 2.1 2.2];
Exp.mwFreq = 9.5;
Exp.SampleFrame = [0 0.3 0.1; 0.2 1 0];

Exp.Range = [300 400];
[Ppos,Ipos] = resfields(Sys,Exp);
Exp.Range = [-400 -300];
[Pneg,Ineg] = resfields(Sys,Exp);

ok = areequal(Pneg,-Ppos,1e-10,'abs') && areequal(Ineg,Ipos,1e-10,'rel');
