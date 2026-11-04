function ok = test() %#ok<FNDEF>

% Davies ENDOR (predefined) needs selective microwave pulses (Exp.tp) and
% an inhomogeneous EPR line (Sys.lwAnisoG). Large couplings give peaks at
% nuI -+ A/2; small couplings (A << 1/tp) are suppressed.

Sys.Nucs = '1H';
Sys.A = 6;  % MHz
Sys.g = planck*9.8e9/bmagn/0.35;
Sys.lwAnisoG = [20 20 20];  % MHz
Sys.lwEndor = 0.05;  % MHz

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.Sequence = 'DaviesENDOR';
Exp.tau = 0.3;  % µs
Exp.tp = [0.2 0.1 0.2];  % µs
Exp.Range = [8 22];  % MHz
Exp.nPoints = 351;

Opt.GridSize = 7;
Opt.nOffsets = 41;

[x,y,info] = saffron(Sys,Exp,Opt);
nuI = larmorfrq('1H',Exp.Field);
peaks = x(islocalmax(real(y)) & real(y)>0.3*max(real(y)));

ok(1) = strcmp(info.Method,'subspace');
ok(2) = numel(peaks)==2 && all(abs(peaks-(nuI+[-1 1]*Sys.A/2))<0.1);

% small coupling: suppressed relative to the large one
Sys.A = 0.4;  % MHz
[~,ySmall] = saffron(Sys,Exp,Opt);
ok(3) = max(real(ySmall)) < 0.2*max(real(y));

% ideal pulses are rejected
Exp.tp = 0;
try
  saffron(Sys,Exp,Opt);
  ok(4) = false;
catch
  ok(4) = true;
end
