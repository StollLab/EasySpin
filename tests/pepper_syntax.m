function ok = test()

% Calling syntax

Sys = struct('S',1/2,'g',[2 2 2.2],'lw',3);
Exp = struct('mwFreq',9.5,'Range',[280 350]);
Opt = struct('Verbosity',0);

spc = pepper(Sys,Exp);
[B,spc] = pepper(Sys,Exp);
[B,spc,out] = pepper(Sys,Exp);

ok(1) = isstruct(out);

spc = pepper(Sys,Exp,Opt);
[B,spc] = pepper(Sys,Exp,Opt);
[B,spc,out] = pepper(Sys,Exp,Opt);

spc = pepper(Sys,Exp,[]); % empty options structure

spc = pepper(Sys,Exp,struct); % empty options structure

ok(2) = true;
