function [err,data] = test(opt,olddata)

% Check whether resfreqs_matrix handles AStrain correctly.

clear Sys Exp
Sys.g = 2 + rand(1,3)*0.1;
Sys.lw = 10;
Sys.A = rand(1,3)*200;
Sys.Nucs = '1H';

Exp.Field = 350;
Opt.PerturbOrder = 'matrix';
Opt.Threshold = 0.001;
[P0,I0] = resfreqs_matrix(Sys,Exp,Opt);

Opt.PerturbOrder = 1;
[P1,I1] = resfreqs_perturb(Sys,Exp,Opt);
err1 = any([sort(P0)-sort(P1)]/1e3>1e-3);

Opt.PerturbOrder = 2;
[P2,I2] = resfreqs_perturb(Sys,Exp,Opt);
err2 = any([sort(P0)-sort(P2)]/1e3>1e-3);

err = err1 || err2;

data = [];
