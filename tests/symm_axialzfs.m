function [err,data] = test(opt,olddata)

%======================================================
% S=1, axial zero-field splitting
%======================================================
Sys = struct('S',1,'g',2,'D',[100 0]);
for n=1:2
  G = symm(Sys);
  err(n) = ~strcmp(G,'Dinfh');
  Sys.ZB11.l = 0; Sys.ZB11.vals = 1;
end
data = [];
