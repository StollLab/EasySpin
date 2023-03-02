% single-crystal spectra with sample rotation
%================================================================

clear, clf

% Substitutional nitrogen center (P1) in diamond
P1.g = 2.0024;
P1.Nucs = '14N';
P1.A = [81 114];             % MHz
P1.lwpp = 0.03;              % mT

% Experimental parameters
Exp.mwFreq = 9.5;            % GHz
Exp.CenterSweep = [339 10];  % mT

% P1 and crystal orientation, crystal symmetry
ma = 54.7361;                               % magic angle (deg)
Exp.MolFrame = [45 ma 0]*pi/180;            % mol. frame ori. of P1 in crystal
Exp.SampleFrame = [0 -ma 0]*pi/180;         % crystal ori. in spectrometer
Exp.CrystalSymmetry = 'Fd-3m';              % space group of diamond (#227)

% Sample rotation axis and angle
nRot = 'x';                  % = x-axis of lab frame (xL)
rho = deg2rad(0:10:180);     % rotate in steps over 180 degrees (half turn)

for k = 1:numel(rho)
  Exp.SampleRotation = {rho(k),nRot};
  [B,spc(k,:)] = pepper(P1,Exp);
end

stackplot(B,spc,3);
xlabel('magnetic field (mT)');
