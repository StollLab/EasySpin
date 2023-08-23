function ok = test()

% Syntax test

Sys = struct('g',[2.01 2.02 2.03],'Nucs','14N','A',[10 20 30]);
Sys.tcorr = 1e-8; Sys.lw = 0.1;
Exp = struct('mwFreq',9.7);
Opt = struct;

spc = chili(Sys,Exp);
[B,spc] = chili(Sys,Exp);
[B,spc,out] = chili(Sys,Exp);

ok(1) = all(size(B)==size(spc));
ok(2) = isstruct(out);

spc = chili(Sys,Exp,Opt);
[B,spc] = chili(Sys,Exp,Opt);
[B,spc,out] = chili(Sys,Exp,Opt);

ok(3) = all(size(B)==size(spc));
ok(4) = isstruct(out);
