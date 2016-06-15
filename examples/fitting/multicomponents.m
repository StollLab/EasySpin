% Fitting a multi-component spectrum
%========================================================================

clear
Sys1.g = [2 2.2];
Sys1.lwpp = 1;
Sys1.weight = 1;
Sys2.g = [2.1 2.15];
Sys2.lwpp = 1;
Sys2.weight = 0.3;

Exp.mwFreq = 9.5;
Exp.Range = [300 350];

y1 = pepper(Sys1,Exp);
y2 = pepper(Sys2,Exp);
y = y1 + y2;
y = addnoise(y,40);

Vary1.g = [0.1 0.1];
Vary2.g = [0.1 0.1];

esfit('pepper',y,{Sys1,Sys2},{Vary1,Vary2},Exp);
