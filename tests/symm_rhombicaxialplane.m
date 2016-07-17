function [err,data] = test(opt,olddata)

%======================================================
% Symmetry detection test:
%  one axial and one rhombic tensor, axial axis in a D2h plane
%======================================================
clear Sys
Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.g = [2 2 3];
Sys.A = [3 2 1];
Sys.gFrame = [rand+pi/10 pi/2 0];
for n=1:2
  G = symm(Sys);
  err(n) = ~strcmp(G,'C2h');
  Sys.ZB11.l = 0; Sys.ZB11.vals = 1;
end
data = [];
