function ok = test(opt)

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
yG = yG*(10.^(p.RCAG/20));
yP = yP*sqrt(p.MWPW*1000);
yc = yc*(p.SPTP*1000);
yT = yT/p.STMP;

% All spectra must be identical
thr = 1e-8;
ok = [areequal(y,yc,thr,'rel') areequal(y,yG,thr,'rel') areequal(y,yP,thr,'rel') areequal(y,yT,thr,'rel')];

if opt.Display
  plot(x,y,x,yc,x,yG,x,yP);
  legend('data','conversion time','gain','power');
  legend boxoff
end
