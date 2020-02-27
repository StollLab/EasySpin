function ok = test()

rng(5);

Sys.S = randi(6)/2;
Sys.g = rand(1,3)*4;
Sys.D = rand(1,2)* 10*30e3;

Exp.Temperature = rand(randi(50),1)*300;
Exp.Field = rand(randi(50),1)*1e4;

Opt.Method = 'partitionfunction';
mu_e = curry(Sys,Exp,Opt);

Opt.Method = 'operator';
mu_o = curry(Sys,Exp,Opt);

Exp.Field = 0;
Opt.Method = 'partitionfunction';
Opt.Output = 'chimolT';
chiTe = curry(Sys,Exp,Opt);

Opt.Method = 'operator';
chiTo = curry(Sys,Exp,Opt);

ok = areequal(mu_e,mu_o,1e-4,'abs') && areequal(chiTe,chiTo,1e4,'abs');
