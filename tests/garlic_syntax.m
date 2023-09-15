function ok = test()

Sys = struct('g',2,'Nucs','1H','lw',0.01,'A',10);
Exp = struct('mwFreq',9.7);
Opt = struct();

y = garlic(Sys,Exp);
[x,y] = garlic(Sys,Exp);
[x,y,info] = garlic(Sys,Exp);

ok(1) = isstruct(info);
ok(2) = all(size(x)==size(y));

y = garlic(Sys,Exp,Opt);
[x,y] = garlic(Sys,Exp,Opt);
[x,y,info] = garlic(Sys,Exp,Opt);

ok(3) = isstruct(info);
ok(4) = all(size(x)==size(y));
ok(5) = isfield(info,'resfields');
