% single-crystal spectra with sample rotation
%================================================================

clear, clf

% Substitutional nitrogen center (P1) in diamond
P1.g = 2.0024;
P1.Nucs = '14N';
P1.A = [81 114];             % MHz
P1.lwpp = 0.03;              % mT

% Spectrometer settings
Exp.mwFreq = 9.5;            % GHz
Exp.CenterSweep = [339 10];  % mT

% P1 orientation, crystal orientation, crystal spacegroup
ma = 54.7356;                       % magic angle (deg), = atan(sqrt(2)) in rad
Exp.MolFrame = [45 ma 0]*pi/180;    % P1 molecular frame orientation in crystal (rad)
Exp.CrystalSymmetry = 'Fd-3m';      % space group of diamond (#227)
Exp.SampleFrame = [0 ma 0]*pi/180;  % crystal orientation in spectrometer (rad)

% Sample rotation axis and angle
nRot = 'x';                  % rotate around lab frame x axis (xL)
rho = deg2rad(0:10:180);     % rotate in 10 degree steps over 180 degrees (rad)
Exp.SampleRotation = {nRot,rho};

Opt.separate = 'orientations';
[B,spc] = pepper(P1,Exp,Opt);

labels = compose('%1.0f°',rho*180/pi);
stackplot(B,spc,'none',1,labels);
xlabel('magnetic field (mT)');
