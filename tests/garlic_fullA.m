function ok = test()

Afull = [6 3 2; 2 5 2; 4 3 10]*10;
Aeig = eig(Afull).';
Aiso = mean(Aeig);

Sys.Nucs = '1H';
Sys.lwpp = 0.1;

Exp.mwFreq = 9.5;
Exp.Range = [337 341];

Sys.A = Afull;
[x1,y1] = garlic(Sys,Exp);
Sys.A = Aeig;
[x2,y2] = garlic(Sys,Exp);
Sys.A = Aiso;
[x3,y3] = garlic(Sys,Exp);

ok = areequal(y1,y2) && areequal(y2,y3);


