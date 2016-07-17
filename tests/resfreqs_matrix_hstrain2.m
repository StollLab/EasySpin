function [err,data] = test(opt,olddata)

% Check whether resfreqs_matrix handles HStrain correctly.

Sys.S = 1/2;
Sys.g = [2 2.05 2.1];
Sys.HStrain = [10 40 100];

Exp.Field = 350;
Exp.Range = [9.5 10.5];

for k=1:5
  ang = rand(1,2)*2*pi;
  z = ang2vec(ang(1),ang(2));
  Exp.CrystalOrientation = [ang(1) ang(2) 0];
  [dum,dum,Wdat] = resfreqs_matrix(Sys,Exp);
  Wdat0 = sqrt(z(:).^2.'*Sys.HStrain(:).^2);
  err(k) = ~areequal(Wdat,Wdat0,1e-4);
end

err = any(err);
data = [];
