function ok = test()

% S=1/2, I=1/2

clear
Sys.Nucs='13C';
Exp.Field=352.3;
Exp.mwFreq=9.5;
Exp.SampleFrame = [0 pi/4 0];

% weak coupling
Sys.A = [5 9];
a=endorfrq(Sys,Exp);
b=endorfrq_perturb(Sys,Exp);
ok1 = all(abs(sort(a)-sort(b))<=0.01);

% strong coupling
Sys.A = [30 40];
a=endorfrq(Sys,Exp);
b=endorfrq_perturb(Sys,Exp);
ok2 = all(abs(sort(a)-sort(b))<=0.1);

ok = ok1 && ok2;

