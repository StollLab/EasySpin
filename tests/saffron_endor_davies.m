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

lastwarn('');
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

% no warning in the regular runs
ok(5) = isempty(lastwarn);

% without Exp.mwFreq and Sys.lwAnisoG, the inversion pulse is centered on
% the EPR line and inverts both hyperfine lines of a spin-1/2 nucleus
% equally: flat spectrum, with a warning
SysC = rmfield(Sys,'lwAnisoG');
SysC.A = 6;  % MHz
ExpC = rmfield(Exp,'mwFreq');
lastwarn('');
evalc('[~,yFlat] = saffron(SysC,ExpC,Opt);');
[~,warnId] = lastwarn;
ok(6) = strcmp(warnId,'EasySpin:saffron:ENDORFlat') && ...
  max(abs(yFlat-yFlat(1)))<1e-10*abs(yFlat(1));

% 14N (I = 1): the lines at mI = 0 and +-1 are inverted differently, so the
% spectrum is not flat, and there is no warning
SysC.Nucs = '14N';
ExpC.Range = [0 8];  % MHz
lastwarn('');
[~,yN] = saffron(SysC,ExpC,Opt);
ok(7) = isempty(lastwarn) && max(abs(yN-yN(1)))>0.1*abs(yN(1));

% ideal pulses are rejected
Exp.tp = 0;
try
  saffron(Sys,Exp,Opt);
  ok(8) = false;
catch
  ok(8) = true;
end
