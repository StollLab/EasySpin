function ok = test()

% Check whether resfreqs_matrix handles AStrain correctly.

clear Sys Exp
Sys.S = 1/2;
Sys.g = 2;
Sys.Nucs = '63Cu';
Sys.A = 100;
Sys.AStrain = 10;

Exp.Field = 350;

Opt.Threshold = 1e-3;

[dum,dum2,Wdat] = resfreqs_matrix(Sys,Exp,Opt);

I = nucspin(Sys.Nucs);
mI = I:-1:-I;
Wdat0 = Sys.AStrain*abs(mI(:));

ok = areequal(Wdat,Wdat0,1e-2,'abs');
