function ok = test()

thisFileName = 'orca/dioxygen_gD.oof';

Sys = orca2easyspin(thisFileName);

D0 = [-0.250117, -0.250117, 1.858352]; % cm^-1
D = Sys.D*1e6/clight/100; % MHz -> cm^-1

ok = areequal(D,D0,1e-4,'abs');
