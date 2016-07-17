function [err,data] = test(opt,olddata)

% Test scaling functionality for BES3T files (.DTA/.DSC)

file = 'eprfiles/E580_Xepr26b6_cwX_WillMyers.DTA';

% Load data, scaled and unscaled
yG = eprload(file,'G'); % receiver gain
yP = eprload(file,'P'); % microwave power
yc = eprload(file,'c'); % conversion/sampling time
yT = eprload(file,'T'); % temperature
%yn = eprload(file,'n'); % number of scans
[x,y,p] = eprload(file);

% Undo the scaling
yG = yG*(10.^(p.RCAG/10));
yP = yP*sqrt(p.MWPW*1000);
yc = yc*(p.SPTP*1000);
yT = yT/p.STMP;

% All spectra must be identical
e = max(abs(y))*1e-8;
ok = [areequal(y,yc,e), areequal(y,yG,e) areequal(y,yP,e) areequal(y,yT,e)];
err = any(~ok);

if opt.Display
  plot(x,y,x,yc,x,yG,x,yP);
end

data = [];
