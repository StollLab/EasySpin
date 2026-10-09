function ok = test()

% Second-order perturbation theory warns for integer S with nuclei (degenerate mS=0 nuclear sublevels)

Sys.S = 1;
Sys.Nucs = '14N';
Sys.A = [10 10 40];
Sys.D = 500;
Exp.Field = 3400;
Opt.PerturbOrder = 2;

lastwarn('');
evalc("resfreqs_perturb(Sys,Exp,Opt);");
ok(1) = ~isempty(lastwarn);

lastwarn('');
Opt.PerturbOrder = 1;
evalc("resfreqs_perturb(Sys,Exp,Opt);");
ok(2) = isempty(lastwarn);

lastwarn('');
Sys.S = 3/2;
Opt.PerturbOrder = 2;
evalc("resfreqs_perturb(Sys,Exp,Opt);");
ok(3) = isempty(lastwarn);
