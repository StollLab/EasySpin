function ok = test()

%======================================================
% D vs. Stevens operators
%======================================================
% Was a bug introduced in 2.2.0 and removed in 2.5.0
% Reported by Frank Loncke 16-Feb-2006

Sys.S = 1;
Sys.g = [1,1,1]*2;
b20 = 100;
b22 = 150;
Sys.B2 = [b22 0 b20 0 0];

Sys1.S = Sys.S;
Sys1.g = Sys.g;
Sys1.D = [-1,-1,2]*b20 + [+1,-1,0]*b22;
S1 = hamsymm(Sys);
S2 = hamsymm(Sys1);
ok(1) = strcmp(S1,'D2h') && strcmp(S2,'D2h');


Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
Sys1.ZB11 = Sys.ZB11;
S1 = hamsymm(Sys);
S2 = hamsymm(Sys1);
ok(2) = strcmp(S1,'D2h') && strcmp(S2,'D2h');
