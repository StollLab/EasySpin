function ok = test()

% S=1, axial zero-field splitting

Sys = struct('S',1,'g',2,'D',[100 0]);
ok(1) = strcmp(symm(Sys),'Dinfh');

Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
ok(2) = strcmp(symm(Sys),'Dinfh');
