% Ground-state triplet of triphenylbenzene di-anion
%----------------------------------------------------------
% see
% - Weil/Bolton, chapter 6.3.3, Fig. 6.6
% - Jesse et al, Mol.Phys. 6, 633 (1963)
%   https://doi.org/10.1080/00268976300100741

clear, clf

Sys.S = 1;
Sys.g = 2;     % isotropic g
Sys.lwpp = 4;  % Gaussian FWHM, mT

% Zero-field splitting in terms of D and E
D = 0.0435;  % cm^-1
E = 0;   % cm^-1
Sys.D = unitconvert([D E],'cm^-1->MHz'); % cm^-1 -> MHz

Exp.mwFreq = 9.150;     % GHz
Exp.Range = [100 500];  % mT

pepper(Sys,Exp); 
