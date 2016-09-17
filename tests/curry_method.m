function [err,data] = test(opt,olddata)

Sys.S =randi(6)/2;
Sys.g = rand(1,3)*4;
Sys.D = rand(1,2)* 10*30e3;

Exp.Temperature = rand(randi(50),1)*300;
Exp.Field = rand(randi(50),1)*1e4;

Opt.Method = 'energies';
me = curry(Sys,Exp,Opt);

Opt.Method = 'operator';
mo = curry(Sys,Exp,Opt);

Exp.Field = 0;
Opt.Method = 'energies';
Opt.Output = 'chiTCGS';
chite = curry(Sys,Exp,Opt);

Opt.Method = 'operator';
chito = curry(Sys,Exp,Opt);

err = ~(areequal(me,mo,1e-4) && areequal(chite,chito,1e4));
data =[];

