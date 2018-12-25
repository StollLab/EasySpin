% trajectories of stochastic rotational diffusion
%==========================================================================
% This script simulates rotational diffusion using the provided rotational
% correlation time. By calling stochtraj_diffusion without any output 
% arguments, rotational trajectories of the body-fixed X-, Y-, and Z-axes 
% are plotted on the unit sphere.

clear

Sys.tcorr = 1e-9;              % rotational correlation time, s

Par.OriStart = [0;0;0];        % starting orientation on the north pole

stochtraj_diffusion(Sys,Par)