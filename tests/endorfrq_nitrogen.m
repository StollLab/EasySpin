function [err,data] = test(opt,olddata)

%-------------------------------------------------------------
% Test 3: Check ENDOR frequencies for 14N system
%-------------------------------------------------------------
Sys = struct('S',1/2,'g',[2 2 2],...
  'Nucs','14N',...
  'A',[3 4 5],...
  'Q',[-1 2 -1],'QFrame',[0 -pi/4 0]);
Par = struct('Field',350,'CrystalOrientation',[0 0 0]);

p0 = [5.378889503898790; 2.379528042419224; 3.567627084365086; 4.134616989405913; 0.566989905040828; 7.758417546318015];
p = endorfrq(Sys,Par);
err = ~areequal(sort(p),sort(p0),1e-7,'rel');
data = [];
