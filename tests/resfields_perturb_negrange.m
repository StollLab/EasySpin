function ok = test()

% Resonance fields at negative fields are the mirror images of those at
% positive fields (perturbation theory)

Sys.g = [2 2.1 2.2];
Sys.Nucs = '1H';
Sys.A = [10 20 30];
Exp.mwFreq = 9.5;
Exp.SampleFrame = [0 0.3 0.1];

Exp.Range = [300 400];
[Ppos,Ipos] = resfields_perturb(Sys,Exp);
Exp.Range = [-400 -300];
[Pneg,Ineg] = resfields_perturb(Sys,Exp);
Exp.Range = [-400 400];
[Pboth,Iboth] = resfields_perturb(Sys,Exp);

ok = areequal(Pneg,-Ppos,1e-10,'abs') && areequal(Ineg,Ipos,1e-10,'rel') && ...
  areequal(Pboth,[Ppos;-Ppos],1e-10,'abs') && areequal(Iboth,[Ipos;Ipos],1e-10,'rel');
