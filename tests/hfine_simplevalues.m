function ok = test()

% explicit comparison

Sys = struct('S',[1/2],'Nucs','1H','g',[2 2 3],'A',[10 12 -44]);

H0 = [-11 0 0 -0.5; 0 11 5.5 0; 0 5.5 11 0; -0.5 0 0 -11];
H1 = hfine(Sys);

ok = areequal(H1,H0,1e-10,'rel');
