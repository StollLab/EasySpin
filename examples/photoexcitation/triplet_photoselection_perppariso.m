% EPR spectrum of triplet photoexcited with linearly polarized light
%===============================================================================
% Example showing that the weighted sum of EPR spectra resulting from
% photoexcitation with E-field perpendicular and parallel to the static 
% magnetic field gives the spectrum expected for complete isotropic 
% light excitation.
%
% Spin system parameters from:
%  - Redman et al., J. Phys. Chem. C, 125, 11782 (2021).
%    https://doi.org/10.1021/acs.jpcc.1c03278
% 

clear; clc; clf;

% Spin system: triplet state acceptor porphyrin
Sys.S = 1;
Sys.g = [2.0041 2.0038 2.0006];
Sys.D = [1033 -272];  % MHz
Sys.lwpp = 2.5;  % mT
Sys.initState = {[0.10 0 0.90],'xyz'};

% Optical transition dipole moment
Sys.tdm = [-3 74]*pi/180;

% Experimental parameters
Exp.mwFreq = 9.5;  % GHz
Exp.Range = [280 400];  % mT
Exp.Harmonic = 0;

% Calculate spectrum for isotropic excitation
Exp.lightBeam = '';
[~,spc_iso] = pepper(Sys,Exp);

% Calculate spectra for light polarized perpendicular and parallel to B0
Exp.lightBeam = 'perpendicular';
[~,spc_perp] = pepper(Sys,Exp);
Exp.lightBeam = 'parallel';
[B,spc_para] = pepper(Sys,Exp);

% Weighted sum of spc_para and spc_perp to give isotropic spectrum
spc_sum = spc_para + 2*spc_perp;

% Plot
hold on; box on;
plot(B,spc_perp)
plot(B,spc_para)
plot(B,spc_sum,'k')
plot(B,spc_iso,'--','Color',[1 1 1]*0.8)
xlabel('{\itB} (mT)')
legend('$\perp$','$\parallel$','$\parallel + 2\perp$','isotropic','interpreter','latex')
