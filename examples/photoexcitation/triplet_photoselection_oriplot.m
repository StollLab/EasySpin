% Photoselection for spin-polarized triplet EPR spectrum
%===============================================================================

% This example illustrates the effect of photoselection on a 
% triplet EPR spectrum for excitation at wavelengths corresponding
% to two perpendicular transition dipole moments in the xy plane
% of the molecular frame (e.g. Qx and Qy in a porphyrin)

clear, clc, clf

% Define spin system
Triplet.S = 1;
Triplet.D = [1045 -282];  % MHz
Triplet.lwpp = 3;  % mT
Triplet.initState = {[0.1 0 0.9],'xyz'};  % non-equilibrium populations

% Define experimental parameters
Exp.mwFreq = 9.5;  % GHz
Exp.Range = [280 400];  % mT
Exp.Harmonic = 0;  % for transient EPR, turn off field modulation

% Light excitation direction and polarization angle (alpha)
k = 'y';
alpha = [0 pi/2];

% Calculate triplet spectra for excitation along two perpendicular optical TDMs
%===============================================================================
Triplet.tdm = 'x';  % orientation of tdm in molecular frame
Exp.lightBeam = {k alpha(1)};
[B,spcx_para] = pepper(Triplet,Exp);
Exp.lightBeam = {k alpha(2)};
[B,spcx_perp] = pepper(Triplet,Exp);

Triplet.tdm = 'y';  % orientation of tdm in molecular frame
Exp.lightBeam = {k alpha(1)};
[B,spcy_para] = pepper(Triplet,Exp);
Exp.lightBeam = {k alpha(2)};
[B,spcy_perp] = pepper(Triplet,Exp);

% Calculate photoselection weights for different orientations
%===============================================================================
% Set up orientational grid
Opt.GridSymmetry = 'C1';
Opt.GridSize = 91;
[gridstruct,tri] = sphgrid(Opt.GridSymmetry,Opt.GridSize);
Orientations = [gridstruct.phi; gridstruct.theta].';
vecs = gridstruct.vecs;

% Calculate weights for all orientations
Triplet.tdm = 'x';
weightsx_para = photoselect(Triplet.tdm,Orientations,k,alpha(1));
weightsx_perp = photoselect(Triplet.tdm,Orientations,k,alpha(2));

Triplet.tdm = 'y';
weightsy_para = photoselect(Triplet.tdm,Orientations,k,alpha(1));
weightsy_perp = photoselect(Triplet.tdm,Orientations,k,alpha(2));

% Plotting
%===============================================================================
subplot(2,4,[1 2])
title('TDM || y_M')
hold on; box on;
plot(B,spcy_para,B,spcy_perp);
legend('parallel','perpendicular')
legend boxoff
xlim(B([1 end]))
ylim([-1 1]*8)
xlabel('{\itB} (mT)')

subplot(2,4,[3 4])
title('TDM || x_M')
hold on; box on;
plot(B,spcx_para,B,spcx_perp);
legend('parallel','perpendicular')
legend boxoff
xlim(B([1 end]))
ylim([-1 1]*8)
xlabel('{\itB} (mT)')

% Plot orientations
subplot(2,4,5)
title('parallel'); hold on; grid on;
trisurf(tri.idx,vecs(1,:),vecs(2,:),vecs(3,:),weightsy_para);
  
subplot(2,4,6)
title('perpendicular'); hold on; grid on;
trisurf(tri.idx,vecs(1,:),vecs(2,:),vecs(3,:),weightsy_perp);

subplot(2,4,7)
title('parallel'); hold on; grid on;
trisurf(tri.idx,vecs(1,:),vecs(2,:),vecs(3,:),weightsx_para);
  
subplot(2,4,8)
title('perpendicular'); hold on; grid on;
trisurf(tri.idx,vecs(1,:),vecs(2,:),vecs(3,:),weightsx_perp);

for i = 5:8
  subplot(2,4,i)
  view([130 40]);
  caxis([0 1])
  shading flat
  axis equal
  xlabel('x_{M}');
  ylabel('y_{M}');
  zlabel('z_{M}');
end


