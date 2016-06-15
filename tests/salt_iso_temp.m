function [err,data] = test(opt,olddata)

% Assure that intensities contain temperature
% dependence when running an isotropic system
% (was a bug in salt, Sep 2012)
  
Sys.Nucs = '13C';
Sys.A = 5;
Sys.lwEndor = 0.5;

Exp.Field = 1150;
Exp.Range = [10 20];

Temps = [5 10 20 40 80];

for t = 1:numel(Temps)
  Exp.Temperature = Temps(t);
  [x,y1(:,t)] = salt(Sys,Exp);
end
S = max(y1);

if opt.Display
  plot(x,y1);
end

data = [];
err = ~(all(diff(S)<0));
