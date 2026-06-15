function ok = test()

% Test reading of ORCA 6.1.0 main output file.
% Molecule: organic radical with O and two H nuclei computed.

Sys = orca2easyspin('orca/model_v610.oof');

ok(1) = Sys.S == 0.5;
ok(2) = strcmp(Sys.Nucs, 'O,H,H');

g0 = [2.0021314 2.0046182 2.0104316];
ok(3) = areequal(sort(Sys.g), sort(g0), 1e-5, 'abs');

Aiso0 = [-26.4807 -20.0047 -19.9016];  % MHz
Aiso = mean(Sys.A, 2).';
ok(4) = areequal(Aiso, Aiso0, 1e-3, 'abs');

end
