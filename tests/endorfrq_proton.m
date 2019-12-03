function [err,data] = test(opt,olddata)

%-------------------------------------------------------------
% Check ENDOR frequencies for 1H system
%-------------------------------------------------------------
Sys = struct('S',1/2,'g',[2 2 2],'Nucs','1H','A',[3 7 12]);
Par = struct('mwFreq',9.5,'Field',350);
Par.CrystalOrientation = [0 0 0];

p = endorfrq(Sys,Par);
p = sort(p);

p0 = [8.902652857; 20.902652857];

err = ~areequal(p,p0,1e-6,'rel');
data = [];
