% Correlated g/A strain in a Cu(II) complex
%==========================================================================
% This example simulates a Cu(II) spectrum including a correlated statistical
% distribution of g and A along the perpendicular direction. The
% correlation manifests a variation of the linewidth across the four
% resolved hyperfine lines at low field.

% See, for example,
%  Froncisz & Hyde, J.Chem. Phys. 73, 3123 (1980), https://doi.org/10.1063/1.440548
%  Giugliarelli & Cannistrato, Chem. Phys. 98, 115 (1985), https://doi.org/10.1016/0301-0104(85)80100-2

clear, clc

Cu.g = [2.08 2.36];

Cu.Nucs = 'Cu';
Cu.A = [20 380];  % MHz

% Define g and A strains, along || direction only
Cu.StrainPars = {'g(2)','A(2)'};
Cu.StrainFWHM = [0.020 50];  % A in MHz
Cu.StrainCorr = -1;  % g-A correlation coefficient, between -1 and 1

% Add additional anisotropic broadening (from unresolved hyperfine etc.)
Cu.HStrain = [90 60];

Exp.mwFreq = 9.6;  % GHz
Exp.Range = [260 350];  % mT

pepper(Cu,Exp);
