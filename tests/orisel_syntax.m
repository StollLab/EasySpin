function ok = test()

% Test 1: syntax check

Sys = struct('S',1/2,'g',[2 2.1 2.2]);
Sys.HStrain = [1 1 1]*50;
Exp = struct('mwFreq',9.5,'Field',320);
Exp.ExciteWidth = 100;
Opt = struct('GridSize',31);

w = orisel(Sys,Exp);
[w,tr] = orisel(Sys,Exp);
w = orisel(Sys,Exp,Opt);
[w,tr] = orisel(Sys,Exp,Opt);

ok = true;
