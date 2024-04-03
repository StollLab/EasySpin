function ok = test()

Sys = struct('S',1/2,'g',[2 3 4],'Nucs','63Cu','A',[56 34 -12]);
Exp = struct('mwFreq',9.5,'Range',[0 600]);
Exp.SampleFrame = rand(1,3)*pi;

p = resfields_eig(Sys,Exp);
[p,i] = resfields_eig(Sys,Exp);
Opt = [];
p = resfields_eig(Sys,Exp,Opt);
Opt = struct('Verbosity',0);
p = resfields_eig(Sys,Exp,Opt);

ok = true;
