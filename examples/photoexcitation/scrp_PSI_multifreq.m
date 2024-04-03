% Spin Correlated Radical Pair 
%===========================================================
% Singlet-born spin-correlated radical pair P700•+ A1•- 
% of Photosystem I at different microwave frequencies
% 
% see
%  - van der Est et al, J. Phys. Chem. B 101, 1437–1443 (1997)
%    https://doi.org/10.1021/jp9622086
%    (Fig. 3, parameters Table 1)
% see also
%  - Poluektov et al, J. Am. Chem. Soc. 127, 11910-11911 (2005)

clear; clc; clf;

% P700•+ A1•- of Photosystem I
Sys.S = [1/2 1/2];
Sys.g = [2.0030 2.0026 2.0023; ... % P700+
         2.0062 2.0051 2.0022];    % A1-
Sys.gFrame = [-10 -128 -83; ...
                0    0   0]*pi/180;

Sys.J = -unitconvert(1e-3,'mT->MHz'); % MHz
Sys.dip = unitconvert(+0.177,'mT->MHz'); % MHz
Sys.eeFrame = [0 90 0]*pi/180;

Sys.lwpp = 0.7; % mT

Sys.initState = 'singlet';

% W-band
Exp(1).Range = [3380 3393.5]; % mT
Exp(1).mwFreq = 95; % GHz

% K-band
Exp(2).Range = [856 870]; % mT
Exp(2).mwFreq = 24.269; % GHz

% X-band
Exp(3).Range = [336.5 350]; % mT
Exp(3).mwFreq = 9.706; % GHz

% Simulate spectra at different microwave frequencies
for i = 1:numel(Exp)  
  Exp(i).nPoints = 4096;
  Exp(i).Harmonic = 0;
  [B{i},spc{i}] = pepper(Sys,Exp(i));
end

% Plotting
for i = 1:numel(Exp)
  subplot(3,1,i)
  plot(B{i},spc{i})
  title(['\nu_{mw} = ',num2str(Exp(i).mwFreq,'%1.2f'),'GHz'])
  xlim(Exp(i).Range)
  xlabel('{\itB} (mT)')
end
