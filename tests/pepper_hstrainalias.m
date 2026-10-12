function ok = test()

% Sys.HStrain (older name) gives the same spectrum as Sys.lwAnisoG

Exp = struct('mwFreq',9.5,'Range',[300 350]);

Sys = struct('g',[2.0 2.05 2.1],'lwAnisoG',[50 80 100]);
[~,y1] = pepper(Sys,Exp);
Sys = struct('g',[2.0 2.05 2.1],'HStrain',[50 80 100]);
warnState = warning('off','EasySpin:HStrainRenamed');
[~,y2] = pepper(Sys,Exp);
warning(warnState);

ok = isequal(y1,y2);
