% Comparison between fast- and slow-motion simulation algorithms
%===============================================================================
% This example shows how diffusion around particular axes of a spin system
% affects the EPR spectrum. For this, we use a rhmobic rotational diffusion
% tensor in Sys.tcorr, with fast diffusion around one axis and very slow
% diffusion around the two others.

clear, clf, clc

% Parameters
%-------------------------------------------------------------------------------
Sys.g = [2 2.05 2.15]; % simple rhombic g tensor, to demonstrate the effects

Exp.mwFreq = 9.5;    % GHz
Exp.Range = [300 360];   % mT
Exp.Harmonic = 0;  % show absorption spectrum

% Here we set the basis size. It needs to be adjusted if spin system parameters
% or diffusion rates are changed. Essentially, increase each number until no
% further changes are seen in the spectra. The numbers must be the higher the
% slower the diffusion.
Opt.LLMK = [20 19 10 10];


% Simulation
%-------------------------------------------------------------------------------
% motional spectra with different diffusion tensors
Sys.tcorr = [0.2 10 10]*1e-9; % fast around x, very slow around y and z
[B,spc1] = chili(Sys,Exp,Opt);

Sys.tcorr = [10 0.2 10]*1e-9; % fast around y, very slow around x and z
[B,spc2] = chili(Sys,Exp,Opt);

Sys.tcorr = [10 10 0.2]*1e-9; % fast around z, very slow around x and y
[B,spc3] = chili(Sys,Exp,Opt);

% ridig-limit spectrum, for reference
[B,spc0] = pepper(Sys,Exp);


% Plotting
%-------------------------------------------------------------------------------
% Plot all simulated spectra. Try to rationalize their shape.
plot(B,spc0,B,spc1,B,spc2,B,spc3);
legend('rigid limit','fast diffusion around x','fast diffusion around y','fast diffusion around z');
xlabel('magnetic field (mT)');
title('Effect of diffusion around individual axes');
