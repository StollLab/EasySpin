% Fitting a multi-component spectrum
%========================================================================

clear
Sys1.g = [2 2.2];
Sys1.lwpp = 1;  % mT
Sys2.g = [2.1 2.15];
Sys2.lwpp = 1;  % mT

Sys1.weight = 1;
Sys2.weight = 0.3;

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [300 350];  % mT

spc1 = pepper(Sys1,Exp);
spc2 = pepper(Sys2,Exp);
spc = spc1 + spc2;
spc = addnoise(spc,40,'n');

Sys1Vary.g = [0.1 0.1];
Sys2Vary.g = [0.1 0.1];

args0 = {{Sys1,Sys2},Exp};
argsvary = {{Sys1Vary,Sys2Vary}};
esfit(spc,@pepper,args0,argsvary);
