% magnetization curve of a cobald dimer complex (crystal)
%=================================================================

clear, clc, clf

% Cobalt dimer with anisotropic g tensors and anisotropic coupling
Co.S = [1/2 1/2];
Co.g = [2 4 7; 2 4 7]; % anisotropic g tensors
Co.gFrame = [0 0 0; 0 0 0]*pi/180;
Co.ee = [3 6 18]*30e3; % coupling, in MHz

% Experiment: temperature and field range
Exp.Temperature = 1.8;
Exp.Field = 0:250:12000;

% Calculate data for powder and three crystal orientations
Exp.CrystalOrientation = []; % powder
magp = curry(Co,Exp);
Exp.CrystalOrientation = [0 0 0]; % crystal, field along crystal z axis
magz = curry(Co,Exp);
Exp.CrystalOrientation = [0 pi/2 0]; % crystal, field along crystal x axis 
magx = curry(Co,Exp);
Exp.CrystalOrientation = [pi/2 pi/2 0]; % crystal, field along crystal y axis
magy = curry(Co,Exp);

% Plotting
B = Exp.Field/1e3;
plot(B,magx,'r',B,magy,'g',B,magz,'b',B,magp,'k');
xlabel('magnetic field (T)');
legend('crystal, B||xC','crystal, B||yC','crystal, B||zC','powder')
legend boxoff
ylabel('longitudinal magnetic moment (\mu_B)');
