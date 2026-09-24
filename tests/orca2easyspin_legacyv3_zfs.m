function ok = test()

% Legacy ORCA 3.0.3 files: zero-field splitting of dioxygen, compared with
% reference values taken manually from dioxygen_gD.oof (cm^-1)

Sys = orca2easyspin('orca/v3.0.3/dioxygen_gD.oof');
D_ref = [-0.250117 -0.250117 1.858352]*100*clight/1e6;  % cm^-1 -> MHz
ok = areequal(Sys.D,D_ref,0.1,'abs');
