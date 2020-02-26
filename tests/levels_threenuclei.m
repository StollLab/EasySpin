function ok = test()

Sys = struct('S',1/2,'g',[2 3 4],'Nucs','14N,14N,63Cu',...
  'A',[4 7 1; 10 23 -8; 78 78 -300]);
phi = 0.4564; theta = 1.14343564;
B = 0;

E = levels(Sys,[phi theta],B);

ok = true;
