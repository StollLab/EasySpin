% Simple isotope mixture
%==========================================================================

clear
Sys.A = [5 50];
Sys.n = [1 1];
Sys.lwpp = 0.05;
Exp.mwFreq = 9.5;
Exp.Range = [330 348];
Exp.nPoints = 1e4;

Sys.Nucs = 'C,11B'; [B,y1] = garlic(Sys,Exp);
Sys.Nucs = 'C,B';  [B,y0] = garlic(Sys,Exp);

plot(B,y1,'g',B,y0,'k');
