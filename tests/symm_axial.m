function [err,data] = test(opt,olddata)

%======================================================
% Test 3: Axial
%======================================================
Sys = struct('S',1/2,'g',[2 2 3]);
for n=1:2
  G = symm(Sys);
  err = ~strcmp(G,'Dinfh');
  Sys.ZB11.l = 0; Sys.ZB11.vals = 1;
end
data = [];
