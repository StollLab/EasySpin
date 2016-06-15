function [err,data] = test(opt,olddata)

%======================================================
% Test 9: Axial+D2h, collinear
%======================================================
Sys = struct('S',1/2,'Nucs','1H','g',[2 2 3],'A',[3 2 1]);

for n=1:2
  G = symm(Sys);
  err(n) = ~strcmp(G,'D2h');
  Sys.ZB11.l = 0; Sys.ZB11.vals = 1;
end
data = [];
