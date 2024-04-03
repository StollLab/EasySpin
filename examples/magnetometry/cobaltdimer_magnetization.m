% magnetization curve of a cobald dimer complex (crystal)
%=================================================================

clear, clc, clf

% Cobalt dimer with anisotropic g tensors and anisotropic coupling
Co.S = [1/2 1/2];
Co.g = [2 4 7; 2 4 7];  % anisotropic g tensors
Co.gFrame = [0 0 0; 0 0 0]*pi/180;
Co.ee = [3 6 18]*30e3;  % coupling, in MHz

% Experiment: temperature and field range
Exp.Temperature = 1.8;  % K
Exp.Field = 0:200:12000;  % mT

% Calculate powder data
magp = curry(Co,Exp);  % powder

% Crystal information
Exp.MolFrame = [0 0 0];  % orientation of dimer within crystal
Exp.CrystalSymmetry = 'P1';  % crystal space group

% Calculate crystal data for three crystal orientations
Exp.SampleFrame = [0 0 0];  % crystal z axis along field
magz = curry(Co,Exp);
Exp.SampleFrame = [0 -pi/2 0];  % crystal x axis along field
magx = curry(Co,Exp);
Exp.SampleFrame = [0 -pi/2 -pi/2];  % crystal y axis along field
magy = curry(Co,Exp);

% Plotting
B = Exp.Field/1e3;
plot(B,magx,'r',B,magy,'g',B,magz,'b',B,magp,'k');
xlabel('magnetic field (T)');
legend('crystal, B0 || x_C','crystal, B0 || y_C','crystal, B0 || z_C','powder')
legend boxoff
ylabel('longitudinal magnetic moment (\mu_B)');
