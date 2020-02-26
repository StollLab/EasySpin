function ok = test()

% Calling syntax

Sys = struct('S',1/2,'g',[2 2 2.2],'lw',3);
Exp = struct('mwFreq',9.5,'Range',[280 350]);
Opt = struct('Verbosity',0);

y = pepper(Sys,Exp);
[x,y] = pepper(Sys,Exp);
[x,y,tr] = pepper(Sys,Exp);
y = pepper(Sys,Exp,Opt);
[x,y] = pepper(Sys,Exp,Opt);
[x,y,tr] = pepper(Sys,Exp,Opt);
x = pepper(Sys,Exp,[]); % empty options structure

ok = true;
