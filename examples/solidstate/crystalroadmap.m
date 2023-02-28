% resonance roadmap from single-crystal rotation
%================================================================

clear, clf

% Spin system parameters
Sys.g = [2, 2.1, 2.2];
Sys.gFrame = [10 20 30]*pi/180;  % rad

% Experimental parameters
Exp.mwFreq = 9.8;  % GHz
Exp.Range = [310 360];  % mT
Exp.CrystalSymmetry = 'P212121';  % space group

% Generate crystal orientations obtained by rotation
rotN = [1 1 0];  % rotation axis
rho = linspace(0,pi,91);  % rotation angle
frame0 = [0 0 0];  % initial crystal orientation
frames = rotateframe(frame0,rotN,rho);
Exp.SampleFrame = frames;

% Calculate resonance fields
Bres = resfields(Sys,Exp);

% plotting
plot(Bres,rho*180/pi);
xlabel('magnetic field (mT)');
ylabel('\rho (deg)');
