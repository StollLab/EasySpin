function [err,data] = test(opt,olddata)

% Assert that HYSCORE sims with product rule give the same intensities are
% sims without product rule.

Sys = struct('Nucs','1H,1H','A',[1 1 3; 4 4 2]);

Exp = struct('Field', 1213.2,...
  'Sequence' , 'HYSCORE',...
  'dt' , 0.004,...
  'nPoints',10,...
  'tau' ,0.1);

Exp.CrystalOrientation = [0 pi/4];

Opt.ProductRule = 0;
[t,y0] = saffron(Sys,Exp,Opt);
Opt.ProductRule = 1;
[t,y1] = saffron(Sys,Exp,Opt);

err = ~areequal(max(abs(y0(:))),max(abs(y0(:))),1e-5);
data = [];
