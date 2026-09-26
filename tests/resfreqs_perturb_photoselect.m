function ok = test()

% Photoselection weights with Sys.tdm

Sys.g = [2 2.1 2.2];
Sys.tdm = 'z';

Exp.Field = 350;
Exp.SampleFrame = [0 0 0];

[~,I0] = resfreqs_perturb(Sys,Exp);
Exp.lightBeam = 'parallel';  % E-field along zL, parallel to tdm
[~,Ipar] = resfreqs_perturb(Sys,Exp);
Exp.lightBeam = 'perpendicular';  % E-field along xL, perpendicular to tdm
[~,Iperp] = resfreqs_perturb(Sys,Exp);
Exp.lightScatter = 0.3;
[~,Iscat] = resfreqs_perturb(Sys,Exp);

ok(1) = Ipar>0;
ok(2) = abs(Iperp)<1e-10*Ipar;
ok(3) = areequal(Iscat,0.3*I0,1e-10,'rel');  % perpendicular weight is zero
