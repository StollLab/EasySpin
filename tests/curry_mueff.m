function [err,data] = test(opt,olddata)

% Make sure the "effective magnetic moment" is correct.
% SI: mueff = sqrt(3*k*T*chi_mol/NA/muB^2/mu0)
% CGS: mueff = sqrt(T*chi*(1/0.12505))

Sys.g = 2;
Sys.S = 1;

Exp.Field = 10;
Exp.Temperature = 100;

Opt.Output = 'mueff';

Opt.Units = 'SI';
mueff_SI = curry(Sys,Exp,Opt);
Opt.Units = 'CGS';
mueff_CGS = curry(Sys,Exp,Opt);

mueff_explicit = Sys.g*sqrt(Sys.S*(Sys.S+1)); % valid only for high-T limit

thr1 = 1e-12;
thr2 = 1e-6;
ok = areequal(mueff_SI,mueff_CGS,thr1,'abs') && ...
     areequal(mueff_SI,mueff_explicit,thr2,'abs');
err = ~ok;

data = [];
