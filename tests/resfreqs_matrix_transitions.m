function ok = test()

% Test transition selection for resfreqs_matrix

clear Sys Exp
Sys.S = 2;
Sys.D = 3*30e3*[1 0];
Sys.lwpp = 2000;

Exp.mwMode = 'perpendicular';
Exp.Temperature = 2;

Opt.Transitions = [1 3];
[nu, int] = resfreqs_matrix(Sys,Exp,Opt);
ok(1) = numel(nu)==1;

Opt.Transitions = [1 3; 2 4];
[nu, int] = resfreqs_matrix(Sys,Exp,Opt);
ok(2) = numel(nu)==2;

Opt.Transitions = [];
[nu, int] = resfreqs_matrix(Sys,Exp,Opt);
ok(3) = numel(nu)==4;

Opt.Transitions = [];
nu = resfreqs_matrix(Sys,Exp,Opt);
ok(4) = numel(nu)==10;
   