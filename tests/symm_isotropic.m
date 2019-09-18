function [err,data] = test(opt,olddata)

%======================================================
% Test 1: Isotropic
%======================================================
Sys = struct('S',1/2,'g',[2 2 2]);
for n=1:2
  G = symm(Sys);
  err(n) = ~strcmp(G,'O3');
  Sys.ZB11.l = 0; Sys.ZB11.vals = 1;
end
data = [];
