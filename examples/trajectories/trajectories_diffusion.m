% trajectories of Markovian jumps between discrete states
%==========================================================================
% This script simulates rotational diffusion using the provided rotational
% correlation time. By calling stochtraj_diffusion without any output 
% arguments, rotational trajectories of the body-fixed X-, Y-, and Z-axes 
% are plotted on the unit sphere.

clear

Sys.tcorr = 1e-9;              % rotational correlation time, s
Sys.Potential = [2,0,0,2.0];   % orientational potential with L=2, M=0, 
                               % K=0, and a coefficient of 2.0

Par.OriStart = [0;0;0];        % starting orientation on the north pole

stochtraj_diffusion(Sys,Par)