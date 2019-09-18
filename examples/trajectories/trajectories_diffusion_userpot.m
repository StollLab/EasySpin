% trajectories of stochastic rotational diffusion
%==========================================================================
% This script simulates rotational diffusion using the provided rotational
% correlation time. By calling stochtraj_diffusion without any output 
% arguments, rotational trajectories of the body-fixed X-, Y-, and Z-axes 
% are plotted on the unit sphere.

clear

Sys.tcorr = 1e-9;              % rotational correlation time, s


N = 70;

abins = N+1;
bbins = N/2+1;
gbins = N+1;

alphaGrid = linspace(0, 2*pi, abins);
betaGrid = linspace(0, pi, bbins);
gammaGrid = linspace(0, 2*pi, gbins);

[AlphaGrid,BetaGrid,GammaGrid] = ndgrid(alphaGrid, betaGrid, gammaGrid);

lambda = 4;

L = 1;
M = 1;
K = 0;
LMK = [L,M,K];
if M==0 && K==0
  Yfun = real(wignerd(LMK,AlphaGrid,BetaGrid,GammaGrid));
else
  Yfun = 2*real(wignerd(LMK,AlphaGrid,BetaGrid,GammaGrid));
end
Sys.Potential = -lambda*Yfun;

Par.OriStart = [0;0;0];        % starting orientation on the north pole
% Par.nTraj = 1;
% Par.nSteps = 10;
Par.dt = Sys.tcorr/20;

stochtraj_diffusion(Sys,Par)