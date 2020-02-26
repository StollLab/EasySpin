function ok = test()

% Test whether salt supports full hyperfine matrices

A = [10 0 0; 0 10 0; 0 0 20];
Sys.A = A;
Sys.Nucs = '1H';
Exp.Field = 350;
Opt.Method = 'perturb2';
[x,y] = salt(Sys,Exp,Opt);
Opt.Method = 'perturb1';
[x,y] = salt(Sys,Exp,Opt);

ok = true;
