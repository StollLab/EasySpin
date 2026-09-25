function ok = test()

% Hybrid method: negative-field resonances including nuclear shifts

Sys.g = [2 2.1 2.2];
Sys.Nucs = '1H';
Sys.A = [10 20 30];
Exp.mwFreq = 9.5;
Exp.SampleFrame = [0 0.3 0.1];
Opt.Method = 'hybrid';

Exp.Range = [300 400];
[Ppos,Ipos] = resfields(Sys,Exp,Opt);
Exp.Range = [-400 -300];
[Pneg,Ineg] = resfields(Sys,Exp,Opt);

ok = areequal(Pneg,-Ppos,1e-10,'abs') && areequal(Ineg,Ipos,1e-10,'rel');
