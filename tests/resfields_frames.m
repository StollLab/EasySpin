function ok = test()

% Compare manually calculated resonance field to resfields() output in the
% presence of non-zero frame tilts (Sys.gFrame, Exp.MolFrame, Exp.SampleFrame)

% This tests 8 of the 9 tilt angles. The test is independent of the 1st angle in
% Exp.SampleFrame.

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [100 400];  % mT

Sys.g = [2 3 4];

Sys.gFrame = [10 20 30]*pi/180;
Exp.MolFrame = [15 25 35]*pi/180;
Exp.SampleFrame = [38 48 67]*pi/180;

Bres = resfields(Sys,Exp);

% Transformation matrices
R_L2S = erot(Exp.SampleFrame);
R_S2M = erot(Exp.MolFrame);
R_M2g = erot(Sys.gFrame);

% Transform field vector from lab frame to g tensor frame
nB0_L = [0; 0; 1];    % B0 is along zLab
nB0_S = R_L2S*nB0_L;
nB0_M = R_S2M*nB0_S;
nB0_g = R_M2g*nB0_M;  % nB0_g is independent of Exp.SampleFrame(1)

% Calculate effective g value and resonance field
g_g = diag(Sys.g);
geff = norm(nB0_g.'*g_g);
Bres_ref = planck*Exp.mwFreq*1e9/bmagn/geff/1e-3;  % mT

ok = all(abs(Bres-Bres_ref)<1e-5);
