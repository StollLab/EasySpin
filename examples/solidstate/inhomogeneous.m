% How to use custom inhomogeneous distributions
%==========================================================================
% In this example, we simulate a simple S=1/2 system with an axial g,
% measured at 9.5 GHz. We assume that the gx value is inhomogeneously
% distributed around a central value, and compute the total broadenend
% spectrum.

clear, clf

% Define the range of gx values and associated weights
%--------------------------------------------------------------------------
N = 100; % number of different systems to simulate and add
gx0 = 2;  % central gx value
gxlw = 0.03;  % linewidth of the distribution
gx = gx0 + linspace(-1,1,N)*gxlw*2.5;  % spread of gx values
weights = gaussian(gx,gx0,gxlw);  % associated probabilities

% System and experiment parameters
%--------------------------------------------------------------------------
Exp.mwFreq = 9.5;  % GHz
Exp.Range = [300 360];  % mT
Sys.lw = 0.8;  % mT
Sys.g = [2, 2.1];

% Simulation loop
%--------------------------------------------------------------------------
spec = 0;      % start with an empty spectrum
for k = 1:N    % run over all possible values of gx
  Sys.g(1) = gx(k);       % set gx value
  [B,thisspec] = pepper(Sys,Exp);   % simulate
  spec = spec + weights(k)*thisspec;  % and add to total, including weight
end

plot(B,spec);
title('Gaussian distribution of gx.');
xlabel('magnetic field (mT)');
