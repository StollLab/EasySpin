function ok = test()

% A powder spectrum must be independent of
% whether D/E or the D diagonal is used to
% provide the zero-field splitting.

D = 99; E = 10;
Sys1.S = 1; Sys2.S = Sys1.S;
Sys1.D = [-1,-1,2]*D/3 + [+1,-1,0]*E;
Sys2.D = [D E];
Sys1.lw = 0.1;
Sys2.lw = Sys1.lw;

Exp.mwFreq = 9.5;
Exp.Range = [333 345];

spc1 = pepper(Sys1,Exp);
spc2 = pepper(Sys2,Exp);

ok = areequal(spc1,spc2,1e-10,'rel');
