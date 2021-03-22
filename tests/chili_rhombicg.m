function ok = test()

% Compare S=1/2 code and general code for rhombic g system

Sys.g = [2.06 2.03 2.0];
Sys.tcorr = 10e-9;

Exp.Field = 333;
Exp.mwRange = [9 10];
Exp.Harmonic = 0;

Opt.LLMK = [6 0 4 4];
Opt.Verbosity = 0;
Opt.MpSymm = false;
Opt.highField = false;
Opt.jKmin = -1;
Opt.evenK = false;

Opt.LiouvMethod = 'fast';
[x1,y1] = chili(Sys,Exp,Opt);

Opt.LiouvMethod = 'general';
[x2,y2] = chili(Sys,Exp,Opt);

y2 = rescaledata(y2,y1,'maxabs');

ok = areequal(y1,y2,1e-2,'abs');
