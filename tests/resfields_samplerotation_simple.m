function ok = test()

Sys.S = 1/2;
Sys.g = [2.0 2.1 2.2];
Exp.mwFreq = 22;
Exp.Range = [0 2000];

Sys.gFrame = [10 -30 70]*pi/180;
Exp.SampleFrame = [40 80 -40];
Exp.MolFrame = [80 10 120]*pi/180;

rotaxis = [1;1;3];
rho = 10*pi/180;
Exp.SampleRotation = {rotaxis,rho};

% Resonance field calculation using resfields
Bres1 = resfields(Sys,Exp);

% Explicit calculation of resonance field
g_g = diag(Sys.g);

R_M2g = erot(Sys.gFrame);
R_g2M = R_M2g.';
g_M = R_g2M*g_g*R_g2M.';

R_S2M = erot(Exp.MolFrame);
R_M2S = R_S2M.';
g_S = R_M2S*g_M*R_M2S.';

R_L2S = erot(Exp.SampleFrame);
R_S2L = R_L2S.';
g_L = R_S2L*g_S*R_S2L.';

Rrot = rotaxi2mat(rotaxis,rho).';  % active rotation!
g_L = Rrot*g_L*Rrot.';

n0_L = [0;0;1];
geff = norm(n0_L.'*g_L);
Bres2 = planck*Exp.mwFreq*1e9/geff/bmagn/1e-3;

ok = areequal(Bres2,Bres1,1e-4,'abs');
