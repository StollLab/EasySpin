function ok = test()

Sys = orca2easyspin('orca\methyl_tilted_all.oof');

% Extract atomic coordinates
xyz = Sys.xyz.';
C = xyz(:,1);
H1 = xyz(:,2);
H2 = xyz(:,3);
H3 = xyz(:,4);

% Calculate molecular z axis
zM = cross(H1-C,H2-C);
zM = zM/norm(zM);

% Assert that z axis of carbon A tensor is along molecular z axis
R_A2M_C = erot(Sys.AFrame(1,:));
zA_C = R_A2M_C(:,3);
ok(1) = areequal(zM,zA_C,1e-3,'rel');

% Assert that z axes of hydrogen A tensors are in molecular xyz plane
R_A2M_H1 = erot(Sys.AFrame(2,:));
zA_H1 = R_A2M_H1(:,3);
R_A2M_H2 = erot(Sys.AFrame(3,:));
zA_H2 = R_A2M_H2(:,3);
R_A2M_H3 = erot(Sys.AFrame(4,:));
zA_H3 = R_A2M_H3(:,3);

ok(2) = areequal(zM.'*zA_H1,0,1e-3,'abs');
ok(3) = areequal(zM.'*zA_H2,0,1e-3,'abs');
ok(4) = areequal(zM.'*zA_H3,0,1e-3,'abs');
