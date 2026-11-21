function ok = test() %#ok<FNDEF>

% Mims ENDOR (predefined), isotropic 1H: the echo is positive and has dips
% at the nuclear transition frequencies nuI -+ A/2, and none at the blind
% spots tau = k/A. The ends of the RF range are far from the lines, so the
% echo there is the echo without RF pulse.

Sys.Nucs = '1H';
Sys.A = 6;  % MHz
Sys.lwEndor = 0.05;  % MHz

Exp.Field = 350;  % mT
Exp.Sequence = 'MimsENDOR';
Exp.tau = 0.1;  % µs
Exp.rfFreq = [8 22];  % MHz
Exp.nPoints = 701;

Opt.GridSize = 7;

[x,y,info] = saffron(Sys,Exp,Opt);
nuI = larmorfrq('1H',Exp.Field);
echoOff = real(y(1));
endor = echoOff - real(y);  % reduction of the echo by the RF pulse
peaks = x(islocalmax(endor) & endor>0.3*max(endor));

ok(1) = strcmp(info.Method,'manifold');
ok(2) = numel(peaks)==2 && all(abs(peaks-(nuI+[-1 1]*Sys.A/2))<0.05);
ok(3) = echoOff>0 && all(endor>-1e-3*max(endor)) && max(endor)>0.05*echoOff;
ok(4) = areequal(y(end),y(1),1e-6,'rel');

% blind spot: flat spectrum, with a warning
Exp.tau = 1/Sys.A;
lastwarn('');
evalc('[~,y] = saffron(Sys,Exp,Opt);');
[~,warnId] = lastwarn;
ok(5) = max(abs(y-y(1)))<1e-10*abs(y(1)) && strcmp(warnId,'EasySpin:saffron:ENDORFlat');
