% spectra from single crystal rotation
%================================================================

clear, clf

% Spin parameters
Sys.g = [2 2.1 2.2];
Sys.gFrame = [10 20 30]*pi/180;
Sys.lwpp = [0.2];

% Experimental parameters
Exp.mwFreq = 9.8;
Exp.Range = [310 360];
Exp.CrystalSymmetry = 'P212121';

% Generate orientations in a single rotation plane
rotN = [1 1 0];  % rotation axis
N = 31;
[phi,theta] = rotplane(rotN,[0 pi],N);
chi = zeros(N,1);
Exp.CrystalOrientation = [phi(:) theta(:) chi];

% Simulate spectra
Opt.Output = 'separate';  % make sure spectra are not added up
[B,spc] = pepper(Sys,Exp,Opt);

% plotting
stackplot(B,spc);
