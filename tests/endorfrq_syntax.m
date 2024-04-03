function ok = test()

Sys = struct('S',1/2,'g',[2 2 2],'Nucs','1H','A',[3 7 12]);
Par = struct('mwFreq',9.5,'Field',350,'SampleFrame',rand(5,3));
Opt = struct('Enhancement','on');

p = endorfrq(Sys,Par);
[p,i] = endorfrq(Sys,Par);
[p,i,t] = endorfrq(Sys,Par);
p = endorfrq(Sys,Par,Opt);
[p,i] = endorfrq(Sys,Par,Opt);
[p,i,t] = endorfrq(Sys,Par,Opt);
[p,i,t,g] = endorfrq(Sys,Par,Opt);

ok = true;
