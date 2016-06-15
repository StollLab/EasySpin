function [err,data] = test(opt,olddata)

%======================================================
% Test 2: Isotropic, tilted frame
%======================================================
Sys = struct('S',1/2,'g',[2 2 2],'gFrame',rand(1,3));
for n=1:2
  G = symm(Sys);
  err = ~strcmp(G,'O3');
  Sys.ZB11.l = 0; Sys.ZB11.vals = 1;
end
data = [];
