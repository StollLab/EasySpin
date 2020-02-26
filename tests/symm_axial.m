function ok = test()

Sys = struct('S',1/2,'g',[2 2 3]);
ok(1) = strcmp(symm(Sys),'Dinfh');
  
Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
ok(2) = strcmp(symm(Sys),'Dinfh');
