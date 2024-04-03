% Comparison of triplet spectrum obtained by photoexcitation with
% unpolarized light with the spectrum expected under isotropic excitation
%===============================================================================
% Example showing that the triplet spectrum recorded with unpolarized light
% in general does not correspond to the spectrum expected for complete
% isotropic excitation (constructed by summation of the spectra obtained
% with light polarized parallel and perpendicular to the field with weights
% of 1 and 2, see also example triplet_photoselection_perppariso.m)
%
% Based on Figure 6 in:
%  - Borovykh et al., J. Phys. Chem. B, 104, 4222 (2000).
%    https://doi.org/10.1021/jp993780a
% 
 
clc, clf, clear

% Spin Hamiltonian parameters and broadening
Sys.S = 1;
Sys.g = [2.002 2.002 2.000];
Sys.D = [563.6 95.9]; % MHz
Sys.lwpp = 1.3; % mT

% T0-populated triplet state
Sys.initState = {[0 1 0],'eigen'};

% Experimental parameters
Exp.mwFreq = 9.5; % GHz
Exp.Range = [310 370]; % mT
Exp.Harmonic = 0;

% Isotropic excitation
[~,spciso] = pepper(Sys,Exp);

% Excitation with unpolarized light
Sys.tdm = [84 82]*pi/180;
Exp.lightBeam = 'unpolarized';
[B,spcunpol] = pepper(Sys,Exp);

% Plot
hold on; box on;
[~,ind] = min(abs(B-321)); % normalize to Z components (as in the paper)
plot(B,spciso/spciso(ind));
plot(B,spcunpol/spcunpol(ind));

legend('isotropically excited','unpolarized','Location','SouthEast')
legend boxoff
xlabel('magnetic field (mT)');
