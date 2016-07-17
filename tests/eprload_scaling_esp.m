function [err,data] = test(opt,olddata)

file = 'eprfiles/100416_wt60min.spc';

% Load data, unscaled and scaled
[x,y,p] = eprload(file);
[x,yn,p] = eprload(file,'n'); % number of scans
[x,yG,p] = eprload(file,'G'); % receiver gain
[x,yP,p] = eprload(file,'P'); % microwave power
[x,yT,p] = eprload(file,'T'); % temperature

% Undo the scaling
yn = yn*p.JSD;
yG = yG*p.RRG;
yP = yP*sqrt(p.MP);
yT = yT/p.TE;

% all spectra must be identical
e = max(abs(y))*1e-8;
ok = [areequal(y,yn,e), areequal(y,yG,e) areequal(y,yP,e) areequal(y,yT,e)];
err = any(~ok);

if opt.Display
  plot(x,y,x,yn,x,yG,x,yP,x,yT);
end

data = [];
