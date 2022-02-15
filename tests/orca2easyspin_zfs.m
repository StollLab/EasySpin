function ok = test()

thisFileName = 'orca/dioxygen_gD.oof';

Sys = orca2easyspin(thisFileName);

D_ref = [-0.250117, -0.250117, 1.858352];  % cm^-1 (as in ORCA output)
D_ref = D_ref*100*clight/1e6;  % cm^-1 -> MHz

ok = areequal(Sys.D,D_ref,1e-4,'abs');
