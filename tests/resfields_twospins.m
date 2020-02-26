function ok = test()

% S=I=1/2 case

Sys = struct('S',1/2,'g',[2 3 4],'Nucs','1H','A',[100 30 60]);
Exp = struct('mwFreq',9.5,'Range',[0 600]);
Exp.CrystalOrientation = [pi/4 pi/4 0];

[p,i] = resfields(Sys,Exp);

p = sort(p);

p0 = [201.7178877473234; 202.1865545783891; 202.5381346934415; 203.0065404415702];
   
ok = all(abs(p-p0)<1e-4);
