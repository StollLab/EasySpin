function [err,data] = test(opt,olddata)
rand('twister',5);

Sys.S =randi(6)/2;
Sys.g = rand(1,3)*4;
Sys.D = rand(1,2)* 10*30e3;

Exp.chiTemperature = rand(randi(50),1)*300;
Exp.chiField = rand(randi(50),1)*1e4;
Exp.mTemperature = rand(randi(50),1)*300;
Exp.mField = rand(randi(50),1)*1e4;

Opt.Output = 'MvsB Chi';

[mC,chiC] = curry(Sys,Exp,Opt);

ExpM.Temperature = Exp.mTemperature;
ExpM.Field = Exp.mField;

Opt.Output = 'MvsB';
mS = curry(Sys,ExpM,Opt);

ExpChi.Temperature = Exp.chiTemperature;
ExpChi.Field = Exp.chiField;

Opt.Output = 'Chi';
chiS = curry(Sys,ExpChi,Opt);

err = ~(areequal(mC,mS) && areequal(chiC,chiS));
data =[];

