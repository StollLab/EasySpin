function ok = test()

% Check interaction between ee/J and hybrid

clear
Sys.S = [0.5,0.5];
Sys.g = [2.00986,2.009967];
Sys.Nucs = ['14N,1H'];
Sys.A = [3.5,0;0,18];
Sys.lwpp = [0 0.05];

Exp.mwFreq = 9.452387;
Exp.Range = [335 337];
Exp.nPoints = 4096;
Exp.SampleFrame = [0 0 0];

Opt.Method = 'hybrid';

Sys.J = 200;
[B,spc1] = pepper(Sys,Exp,Opt);

Sys = rmfield(Sys,'J');
Sys.ee = [200 200 200];
[B,spc2] = pepper(Sys,Exp,Opt);

ok = areequal(spc1,spc2,1e-4,'rel');
