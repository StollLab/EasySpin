function [err,data] = test(opt,olddata)

%======================================================
% Test 8: 2x Axial, tilted general
%======================================================
AFrame = pi*rand(1,3);

Sys = struct('S',1/2,'Nucs','1H','g',[2 2 3],'A',[3 3 10],'AFrame',AFrame);
for n=1:2
  G = symm(Sys);
  % Here the QM method of symm() has difficulties, since it doesn't
  % find the C2h frame. It therefore returns Ci.
  err(n) = ~strcmp(G,'C2h') & ~strcmp(G,'Ci');
  Sys.ZB11.l = 0; Sys.ZB11.vals = 1;
end
data = [];
