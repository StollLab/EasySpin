function ok = test()

S = [1/2 1/2];
g = [2 2.1 2.2];
J = 100;
dip = [-2 1 1]*30;

Sys{1}.S = S;
Sys{1}.g = [g; g];
Sys{1}.ee = J + dip;
Sys{2}.S = S;
Sys{2}.g = [g; g];
Sys{2}.J = J;
Sys{2}.dip = dip;

frame2 = [-15 90 0]*pi/180;
B0 = 400;
R = erot(frame2);
B = B0*[0;0;1];

for k = 1:2
  SysA = Sys{k};
  SysB = Sys{k};
  SysB.gFrame = [frame2; frame2];
  SysB.eeFrame = [frame2];
  
  H1 = ham(SysA,B);
  H2 = ham(SysB,R.'*B);
  
  ok(k) = areequal(eig(H1),eig(H2),1e-10,'abs');
end
