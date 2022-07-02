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

% Experimental parameters
Exp.mwFreq = 9.5; % GHz
Exp.Range = [310 370]; % mT
Exp.Harmonic = 0;

% Construct recombination triplet spectrum from separate transitions
% (see documentation on spin polarization and the example
%  triplet_recombination.m)
Opt.Output = 'separate'; 
% Vector of populations, lowest to highest energy level
Pop = [0 1 0];      % vector of populations for the three levels
Pola = -diff(Pop);  % vector of polarization across the two transitions

% Isotropic excitation
[~,spcsepiso] = pepper(Sys,Exp,Opt);

% Excitation with unpolarized light
Sys.tdm = [84 82]*pi/180;
Exp.lightBeam = 'unpolarized';
[B,spcsepunpol] = pepper(Sys,Exp,Opt);

% EPR spectra are polarization-weighted sum of individual transitions
spciso = Pola(1)*spcsepiso(1,:) + Pola(2)*spcsepiso(2,:); % polarized spectrum
spcunpol = Pola(1)*spcsepunpol(1,:) + Pola(2)*spcsepunpol(2,:); % polarized spectrum##

% Plot
hold on; box on;
[~,ind] = min(abs(B-321)); % normalize to Z components (as in the paper)
plot(B,spciso/spciso(ind));
plot(B,spcunpol/spcunpol(ind));

legend('isotropically excited','unpolarized','Location','SouthEast')
legend boxoff
xlabel('magnetic field (mT)');
