% Photoselection effects on spin-correlated radical pair spectrum
%============================================================================
% Effects of photoselection on the spin-correlated radical pair spectrum of 
% [P•+ QA•-] in Rhodobacter sphaeroides reaction centers
%
% see
% Proskuryakov, I. I., Klenina, I. B., Borovykh, I. V., Gast, P., Hoff, A. J. 
% Chem. Phys. Lett. (1999), 299, 566–570.
% https://doi.org/10.1016/S0009-2614(98)01299-8
%

clear; clc; clf;

% P•+ QA•-
% Spin system parameters from van der Est, A. et al.
% Chem. Phys. Lett. (1993), 212, 561-568.
% https://doi.org/10.1016/0009-2614(93)85487-9
Sys.S = [1/2 1/2];
Sys.g = [2.0033 2.0025 2.0021; ... % P+
         2.0066 2.0054 2.0022];    % QA-
Sys.gFrame = [-242.2 -64 -13.7; ...
                0 0 0]*pi/180;
Sys.J = 0; % MHz
Sys.dip = unitconvert(0.124/3,'mT->MHz'); % MHz
Sys.eeFrame = [-68.5 -71.6 0]*pi/180;
Sys.lwpp = 0.33; % mT

% Singlet-born radical pair
Sys.initState = 'singlet';

% Experimental parameters
Exp.Range = [328 332.5]; % mT
Exp.mwFreq = 9.26; % GHz
Exp.Harmonic = 0;

% Transition dipole moment
Sys.tdm = 'z';

Exp.lightBeam = 'unpolarized';
[B,spc(:,3)] = pepper(Sys,Exp);

Exp.lightBeam = 'perpendicular';
[B,spc(:,2)] = pepper(Sys,Exp);

Exp.lightBeam = 'parallel';
[B,spc(:,1)] = pepper(Sys,Exp);

% Plot spectra
stackplot(B,spc,'none',1,{'parallel','perpendicular','unpolarized'})
xlim(Exp.Range)
xlabel('{\itB} (mT)')

