function ok = test()

% Field range across zero gives resonances at positive and negative fields

Sys.g = [2 2.1 2.2];
Exp.mwFreq = 9.5;
Exp.SampleFrame = [0 0.3 0.1];

Exp.Range = [0 400];
[P0,I0] = resfields(Sys,Exp);
Exp.Range = [-400 400];
[P,I] = resfields(Sys,Exp);

ok = areequal(P,[P0;-P0],1e-10,'abs') && areequal(I,[I0;I0],1e-10,'rel');
