function ok = test()

% Test main ouput file against property output file.

Sys1 = orca2easyspin('orca/dioxygen_singlepoint.oof');
Sys2 = orca2easyspin('orca/dioxygen_singlepoint_property.txt');

ok(1) = all(cellfun(@(a,b)isequal(a,b),fieldnames(Sys1),fieldnames(Sys2)));
ok(2) = areequal(Sys1.g,Sys2.g,1e-6,'rel');
ok(3) = areequal(erot(Sys1.gFrame(1,:)),erot(Sys2.gFrame(1,:)),1e-6,'rel');
ok(4) = areequal(Sys1.D,Sys2.D,1e-6,'rel');
ok(5) = areequal(erot(Sys1.DFrame(1,:)),erot(Sys2.DFrame(1,:)),1e-6,'rel');
ok(6) = areequal(Sys1.A,Sys2.A,1e-6,'rel');
ok(7) = areequal(erot(Sys1.AFrame(1,:)),erot(Sys2.AFrame(1,:)),1e-6,'rel');
ok(8) = areequal(erot(Sys1.AFrame(2,:)),erot(Sys2.AFrame(2,:)),1e-6,'rel');
