% resonance roadmap from single-crystal rotation
%================================================================

clear, clf

% Spin system parameters
Sys.g = [2, 2.1, 2.2];
Sys.gFrame = [10 20 30]*pi/180;  % rad

% Sample parameters
Exp.CrystalSymmetry = 'P212121';     % space group
Exp.MolFrame = [50 20 110]*pi/180;   % orientation of spin system in crystal

% Initial crystal orientation
Exp.SampleFrame = [40 80 20];   % initial crystal orientation, rad

% Information about crystal orientation/rotation
rotaxis = [1 0 0];  % rotation axis
rho = linspace(0,pi,91);  % rotation angle, rad
Exp.SampleRotation = {rotaxis,rho};

% Field and frequency settings
Exp.mwFreq = 9.8;  % GHz
Exp.Range = [310 360];  % mT

% Calculate resonance fields
Bres = resfields(Sys,Exp);

% Plotting
plot(Bres,rho*180/pi);
xlabel('magnetic field (mT)');
ylabel('rotation angle (deg)');
