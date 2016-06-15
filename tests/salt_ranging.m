function [err,data] = test(opt,olddata)

% Manual and automatic rf ranging for ENDOR perturbation theory

clear Sys Exp Exp1 Exp2 Exp3
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [3 6 10];
Sys.lwEndor = 0.055; % in MHz

Exp.Field = 450;
Exp1 = Exp;
Exp2 = Exp; Exp2.Range = [10 50];
Exp3 = Exp; Exp3.CenterSweep = [30 50];

Opt.nKnots = [12 0];

Opt.Method = 'perturb';
[x,y] = salt(Sys,Exp1,Opt);
[x,y] = salt(Sys,Exp2,Opt);
[x,y] = salt(Sys,Exp3,Opt);

Opt.Method = 'matrix';
[x,y] = salt(Sys,Exp1,Opt);
[x,y] = salt(Sys,Exp2,Opt);
[x,y] = salt(Sys,Exp3,Opt);


err = 0;
data = [];
