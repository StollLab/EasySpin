function ok = test()

% Resonance fields at negative fields are the mirror images of those at
% positive fields (eigenfield method)

Sys.g = [2 2.1 2.2];
Exp.mwFreq = 9.5;
Exp.SampleFrame = [0 0.3 0.1];

Exp.Range = [300 400];
[Ppos,Ipos] = resfields_eig(Sys,Exp);
Exp.Range = [-400 -300];
[Pneg,Ineg] = resfields_eig(Sys,Exp);

ok = areequal(Pneg,-Ppos,1e-8,'abs') && areequal(Ineg,Ipos,1e-8,'rel');
