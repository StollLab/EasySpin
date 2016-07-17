function [err,data] = test(opt,olddata)

%======================================================
% Test 7: 2x Axial, z axis tilt general
%======================================================
Sys = struct('S',1/2,'Nucs','1H','g',[2 2 3],'A',[3 3 1],'AFrame',[0 rand+pi/20 0]);
for n=1:2
  G = symm(Sys);
  err(n) = ~strcmp(G,'C2h');
  Sys.ZB11.l=0;Sys.ZB11.vals = 1;
end
data = [];