function ok = test()

% 2x Axial, tilted general
AFrame = pi*rand(1,3);

Sys = struct('S',1/2,'Nucs','1H','g',[2 2 3],'A',[3 3 10],'AFrame',AFrame);
G = hamsymm(Sys);
% Here the QM method of hamsymm() has difficulties, since it doesn't
% find the C2h frame. It therefore returns Ci.
ok(1) = strcmp(G,'C2h') || strcmp(G,'Ci');

Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
G = hamsymm(Sys);
ok(2) = strcmp(G,'C2h') || strcmp(G,'Ci');
