% Simple isotope mixture
%==========================================================================

clear
Sys.A = [5 50];
Sys.n = [1 1];
Sys.lwpp = 0.05;
Exp.mwFreq = 9.5;
Exp.Range = [330 348];
Exp.nPoints = 1e4;

Sys.Nucs = 'C,11B'; [x,y1] = garlic(Sys,Exp);
Sys.Nucs = 'C,B';  [x,y0] = garlic(Sys,Exp);

plot(x,y1,'g',x,y0,'k');
