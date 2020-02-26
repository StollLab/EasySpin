function ok = test()

% calling syntax

Sys = struct('S',1/2,'g',[2 3 4],'Nucs','63Cu','A',[56 34 -12]);
Exp = struct('mwFreq',9.5,'Range',[0 600]);
Ori = [0.1454325 0.75454345 0];
Exp.CrystalOrientation = Ori;

p = resfields(Sys,Exp);
[p,i] = resfields(Sys,Exp);
[p,i,w] = resfields(Sys,Exp);
[p,i,w,tr] = resfields(Sys,Exp);
%[p,i,w,tr,g] = resfields(Sys,Exp);
Opt = [];
p = resfields(Sys,Exp,Opt);
Opt = struct('Verbosity',0);
p = resfields(Sys,Exp,Opt);

ok = true;
