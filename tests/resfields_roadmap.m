function ok = test()

% Check that resfields() returns all sites in a crystal roadmap calculation

clear

% Spin parameters
Sys.g = [2, 2.1, 2.2];
Sys.gFrame = [10 20 30]*pi/180;

% Experimental parameters
Exp.mwFreq = 9.8;
Exp.Range = [310 360];
Exp.CrystalSymmetry = 'P212121';

% Generate orientations in a single rotation plane
rotN_L = [1; 1; 0];         % rotation axis (in lab frame)
rho = linspace(0,pi,31);    % list of rotation angles
frame0 = [0 0 0];           % initial crystal orientation
Exp.SampleFrame = rotateframe(frame0,rotN_L,rho);

% Simulate spectra
Bres = resfields(Sys,Exp);

ok = areequal(size(Bres),[4 numel(rho)]);
