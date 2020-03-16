function ok = test()

% Test whether B2Frame works correctly

% Spin system parameters
S = 2;
D = 1000;
E = D*0.2;

% Define spin systems with D/DFrame and B2/B2Frame
SysD.S = S;
SysD.D = [D E];
SysB.S = S;
SysB.B2 = [E 0 D/3 0 0];

% Test 1: tensor collinear with molecular frame
HD = zfield(SysD);
HB = zfield(SysB);
ok(1) = areequal(HD,HB,1e-10,'rel');

% Test 2: tensor tilted relative to molecular frame, beta-only tilt
ang = [0 pi/7.5 0];
SysB.B2Frame = ang;
SysD.DFrame = ang;
HD = zfield(SysD);
HB = zfield(SysB);
ok(2) = areequal(HD,HB,1e-10,'rel');

% Test 2: tensor tilted relative to molecular frame, general tilt
ang = [pi/3 pi/7.5 pi/5];
SysB.B2Frame = ang;
SysD.DFrame = ang;
HD = zfield(SysD);
HB = zfield(SysB);
ok(3) = areequal(HD,HB,1e-10,'rel');
