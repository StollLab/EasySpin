% CISS effect for a spin-correlated radical pair
% ==========================================================
% Simulation of the influence of chirality-induced spin
% selectivity on the spectrum of a spin-correlated radical
% pair generated through photo-driven charge transfer
% in a DNA hairpin
%
% see
%  - Latawiec, E. I. et al., Proc. Natl. Acad. Sci. 122, e2515120122 (2025)
%    https://doi.org/10.1073/pnas.2515120122
%

clear; clc;

% Spin system
% ==========================================================
Sys. S = [1/2 1/2];
Sys.g = [2.0044 2.00467 2.0022; % NDI•-
         2.0033 2.0043 2.00235];% Sd•+
Sys.gFrame = [  0 0 0;
              204 0 0]*pi/180;
Sys.lwpp = 0.7; % mT

% Spin-spin interactions
r = 2.04; % nm
dip = (mu0/4/pi)*bmagn^2*gfree^2/planck/1e6./(r*1e-9).^3; % MHz
J = 0.005; % MHz

% molecular frame defined as charge separation/chiral axis 
% (= dipolar axis, eeFrame = [0 0 0])
Sys.ee = J + dip*[1 1 -2];

% Experimental settings
% ==========================================================
Exp.mwFreq = 34.75; % GHz
Exp.Range = [1235 1243.5]; % mT
Exp.Harmonic = 0; % direct detection


% Simulation of a singlet-born SCRP
% ==========================================================
Sys.initState = {[0 0 0 1],'coupled'};
[Bsim,simsingletborn] = pepper(Sys,Exp);

clf
hold on; box on;
plot(Bsim,simsingletborn,'k','LineWidth',1,'DisplayName','singlet-born SCRP')

% Spin polarization pattern arising from CISS effect
% ==========================================================
CISSpercentage = [0 25 50 75 100]; % %, = (1-cos(chi))*100

for i = 1:numel(CISSpercentage)

  chi(i) = acos(-(CISSpercentage(i)/100-1));

  % Coherence-based CISS state vector
  pc = [0 1i*sin(chi(i)/2) 0 cos(chi(i)/2)]; % [T+ T0 T- S]
  % Initial state defined through density matrix in coupled basis
  Sys.initState = {pc'*pc,'coupled'};

  [Bsim,simciss(i,:)] = pepper(Sys,Exp);

  plot(Bsim,simciss(i,:),'DisplayName',sprintf('%1.0f%s CISS',CISSpercentage(i),'%'))

end

colororder(parula(numel(CISSpercentage)+2))
legend('Location','SouthEast')
xlim(Exp.Range)
xlabel('magnetic field (mT)')
