function [err,data] = test(opt,olddata)

%======================================================
% Test 12: D2h+D2h, one axis coincides
%======================================================
Sys = struct('S',1/2,'Nucs','1H','g',[2 2.5 3],'A',[3 2 1],'AFrame',[rand+pi/20 0 0]);
for n=0:1
G = symm(Sys);
err(1+n) = ~strcmp(G,'C2h');
Sys.ZB11.l= 0;Sys.ZB11.vals= 1; %isotropic term, no change in symmetry
end
data = [];
