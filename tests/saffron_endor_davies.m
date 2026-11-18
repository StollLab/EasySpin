function ok = test() %#ok<FNDEF>

% Davies ENDOR (predefined) needs selective microwave pulses (Exp.tp) and
% an inhomogeneous EPR line (Sys.lwAnisoG). The echo is negative (inverted),
% and the RF pulse reduces the inversion. Large couplings give peaks at
% nuI -+ A/2; small couplings (A << 1/tp) are suppressed. The ends of the
% RF range are far from the lines, so the echo there is the echo without
% RF pulse.

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

[x,y,info] = saffron(Sys,Exp,Opt);
nuI = larmorfrq('1H',Exp.Field);
endor = real(y-y(1));  % reduction of the inversion by the RF pulse
peaks = x(islocalmax(endor) & endor>0.3*max(endor));

ok(1) = strcmp(info.Method,'subspace');
ok(2) = numel(peaks)==2 && all(abs(peaks-(nuI+[-1 1]*Sys.A/2))<0.1);
ok(3) = real(y(1))<0 && all(endor>-1e-3*max(endor));

% small coupling: suppressed relative to the large one
Sys.A = 0.4;  % MHz
[~,ySmall] = saffron(Sys,Exp,Opt);
ok(4) = max(real(ySmall-ySmall(1))) < 0.2*max(endor);

% ideal pulses are rejected
Exp.tp = 0;
try
  saffron(Sys,Exp,Opt);
  ok(5) = false;
catch
  ok(5) = true;
end
