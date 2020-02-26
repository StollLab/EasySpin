function ok = test(opt,olddata)

% Isotopologues with several equivalent groups of nuclei

a11 = 0.9; a12 = 1 - a11;
a21 = 0.8; a22 = 1 - a21;
Sys.Nucs = '(1,2)H,(12,13)C';
Sys.n = [2 2];
Sys.Abund = {[a11 a12],[a21 a22]};

ab0 = [a11^2*a21^2, a11^2*a21*a22*2, a11^2*a22^2, 2*a11*a12*a21^2, ...
  4*a11*a12*a21*a22, 2*a11*a12*a22^2, a12^2*a21^2, 2*a12^2*a21*a22, a12^2*a22^2];

Iso = isotopologues(Sys,0);

ok = areequal(ab0,[Iso.weight],1e-7,'abs');
