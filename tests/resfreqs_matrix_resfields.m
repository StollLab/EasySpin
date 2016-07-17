function [err,data] = test(opt,olddata)

% Check whether resfreqs_matrix returns the microwave frequency if
% supplied with resonance fields determined by resfreqs.

Sys.S = 1;
Sys.D = 30e3*0.2*[1 0.2];
Sys.Nucs = '1H';
Sys.A = [100 150 -200];

Exp1.mwFreq = 10;
Exp1.Range = [150 400];

for iTest=1:6
  Exp1.CrystalOrientation = rand(1,3)*pi;
  Exp2.CrystalOrientation = Exp1.CrystalOrientation;
  B = resfields(Sys,Exp1);
  for iB = 1:numel(B)
    Exp2.Field = B(iB);
    nu = resfreqs_matrix(Sys,Exp2)/1e3;
  end
  err(iTest) = ~any(abs(nu-Exp1.mwFreq)<1e-6);
end

err = any(err);

data = [];
