function [err,data] = test(opt,olddata)

rand('twister',5);

Sys.S = randi_(6)/2;
Sys.g = rand(1,3)*4;
Sys.D = rand(1,2)* 10*30e3;

Exp.Temperature = rand(randi_(50),1)*300;
Exp.Field = rand(randi_(50),1)*1e4;

Opt.Method = 'energies';
mu_e = curry(Sys,Exp,Opt);

Opt.Method = 'operator';
mu_o = curry(Sys,Exp,Opt);

Exp.Field = 0;
Opt.Method = 'energies';
Opt.Output = 'chimolT';
chiTe = curry(Sys,Exp,Opt);

Opt.Method = 'operator';
chiTo = curry(Sys,Exp,Opt);

err = ~(areequal(mu_e,mu_o,1e-4) && areequal(chiTe,chiTo,1e4));
data =[];
