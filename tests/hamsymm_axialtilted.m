function ok = test()

% Axial, tilted frame

Sys = struct('S',1/2,'g',[2 2 3],'gFrame',rand(1,3));
ok(1) = strcmp(hamsymm(Sys),'Dinfh');

Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
ok(2) = strcmp(hamsymm(Sys),'Dinfh');
