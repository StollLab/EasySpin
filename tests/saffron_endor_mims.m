function ok = test() %#ok<FNDEF>

% Mims ENDOR (predefined), isotropic 1H: peaks at the nuclear transition
% frequencies nuI -+ A/2, and blind spots at tau = k/A.

Sys.Nucs = '1H';
Sys.A = 6;  % MHz
Sys.lwEndor = 0.05;  % MHz

Exp.Field = 350;  % mT
Exp.Sequence = 'MimsENDOR';
Exp.tau = 0.1;  % µs
Exp.Range = [8 22];  % MHz
Exp.nPoints = 701;

Opt.GridSize = 7;

[x,y,info] = saffron(Sys,Exp,Opt);
nuI = larmorfrq('1H',Exp.Field);
peaks = x(islocalmax(real(y)) & real(y)>0.3*max(real(y)));

ok(1) = strcmp(info.Method,'subspace');
ok(2) = numel(peaks)==2 && all(abs(peaks-(nuI+[-1 1]*Sys.A/2))<0.05);
ok(3) = all(real(y)>-1e-3*max(real(y)));  % positive ENDOR effect

% blind spot
Exp.tau = 1/Sys.A;
[~,y] = saffron(Sys,Exp,Opt);
ok(4) = max(abs(y))<1e-12;
