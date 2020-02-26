function ok = test()

% Test Options.Enhancement

Sys = struct('S',1/2,'g',[2 2 2],'Nucs','1H','A',[3 7 12]);
Par = struct('Field',350,'CrystalOrientation',[pi/5 pi/7 pi/6]);
Opt = struct('unused',NaN);

Opt.Enhancement = 'off';
[p0,i0] = endorfrq(Sys,Par,Opt);
Opt.Enhancement = 'on';
[p1,i1] = endorfrq(Sys,Par,Opt);

ok = abs(diff(i0))<=abs(diff(i1));
