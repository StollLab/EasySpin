function [err,data] = test(opt,olddata)

Sys.S = [1/2 1];
Sys.g = [2 3.5];
Sys.D = [0 0; [-1 2]*99];

Sys.J = 50*30e3;

F = spinladder(Sys);

ok = numel(F)==2;
F1 = F{1};
F2 = F{2};
ok = ok && F1.S==0.5 && F2.S==1.5;
ok = ok && F1.g==4 && F2.g==3;
ok = ok && all(F2.D==[77 -55 -22]);

err = ~ok;

data = [];
