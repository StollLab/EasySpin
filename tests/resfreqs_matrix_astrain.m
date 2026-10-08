function ok = test()

% Check whether resfreqs_matrix handles A strain correctly.

clear Sys Exp
Sys.S = 1/2;
Sys.g = 2;
Sys.Nucs = '63Cu';
Sys.A = 100;
Sys.StrainPars = {'A'};
Sys.StrainFWHM = 10;

Exp.Field = 350;

Opt.Threshold = 1e-3;

[dum,dum2,Wdat] = resfreqs_matrix(Sys,Exp,Opt);

I = nucspin(Sys.Nucs);
mI = (I:-1:-I).';
nu = Sys.g*bmagn*Exp.Field*1e-3/planck/1e6; % electron Zeeman frequency, MHz

% Derivative of the transition frequency with respect to A, up to second
% order in A/nu: mI + A/nu*(I(I+1)-mI^2)
Wdat0 = Sys.StrainFWHM*abs(mI + Sys.A/nu*(I*(I+1)-mI.^2));

ok = areequal(Wdat,Wdat0,0.01,'abs');
