% Plot of (D,E) distribution given by D and correlated D/E strains
%===============================================================================
clear, clc

% This script plots the Gaussian distribution over zero-field splitting
% parameter D and E that EasySpin models based on the spin system fields
% Sys.D and the strain fields Sys.StrainPars, Sys.StrainFWHM and Sys.StrainCorr.
% The values are typical of a spin triplet (S = 1) such as an organic triplet
% or biradical.

% Relevant spin system parameters
% (Change these parameters to see how they affect the distribution)
% (The correlation coefficient should be between -1 and 1.)
Sys.D = [1500 150]; % D and E, in MHz
Sys.StrainPars = {'D(1)','D(2)'}; % strains of D and E
Sys.StrainFWHM = [200 50]; % FWHM of D and E Gaussian distributions, MHz
Sys.StrainCorr = -0.9; % correlation coefficient between D and E

% Get center of (D,E) distribution
Dcenter = Sys.D(1);
Ecenter = Sys.D(2);
mu = [Dcenter Ecenter];

% Get widths (fwhm) of (D,E) distribution and convert to standard deviations
Dfwhm = Sys.StrainFWHM(1); % FWHM of Gaussian distribution of D
Efwhm = Sys.StrainFWHM(2); % FWHM of Gaussian distribution of E
Dsigma = Dfwhm/sqrt(2*log(2))/2; % convert FWHM to standard deviation
Esigma = Efwhm/sqrt(2*log(2))/2; % convert FWHM to standard deviation

% Get the correlation coefficient between D and E
rDE = Sys.StrainCorr;

% Build the covariance matrix
R12 = rDE*Dsigma*Esigma; % off-diagonal element of covariance matrix
CovMatrix = [Dsigma^2 R12; R12 Esigma^2];

% Generate appropriate ranges of D and E values
D = Dcenter + Dsigma*linspace(-3,3,501);
E = Ecenter + Esigma*linspace(-3,3,502);

% Generate 2D arrays for D and E, and evaluate the 2D Gaussian over (D,E)
[D_,E_] = meshgrid(D,E);
weights = gaussian2d(D_,E_,mu,CovMatrix);

% Plotting
pcolor(D,E,weights);
hold on
contour(D,E,weights);
set(gca,'Layer','top');
shading flat
grid on
xlabel('D (MHz)');
ylabel('E (MHz)');
colorbar
title(sprintf('Gaussian (D,E) distribution, correlation coefficient %g',Sys.StrainCorr));

% Function that calculates a two-dimensional Gaussian distribution over
% x1 and x2 centered at mu and with covariance matrix Sigma
function y = gaussian2d(x1,x2,mu,Sigma)
y = zeros(size(x1));
x = [x1(:) x2(:)];
dx = x-mu(:).';
invSigma = inv(Sigma);
for i = 1:numel(y)
  y(i) = exp(-dx(i,:)*invSigma*dx(i,:).'/2);
end
y = y/sqrt(det(Sigma))/sqrt((2*pi)^2);
end
