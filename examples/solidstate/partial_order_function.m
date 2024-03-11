% Partial ordering using custom function
%==========================================================================
% This example shows how to provide an orientational distribution function
% via Exp.Ordering to specify partial ordering for a spin center.

clear, clf, clc

% Spin center
Sys.g = [2 2.08 2.2];
Sys.lwpp = 0.8;  % mT

% Experimental parameters
Exp.mwFreq = 9.5;  % GHz
Exp.Range = [300 350];  % mT
Exp.Harmonic = 0;  % absorption spectra, for visual clarity

% Options
Opt = struct;

% (0) Isotropic orientational distribution
[B,spec_iso] = pepper(Sys,Exp,Opt);

% Distribution widths along orientational angles
% (beta and gamma are the second and third of the Euler angles that
% transform the sample frame into the molecular frame)
w_beta = deg2rad(30);
w_gamma = deg2rad(30);

% Elementary distribution function for an angle x, centered at x0, width w
f = @(x,x0,w) exp(-((x-x0)./w).^2);

% (1) Partially oriented, molecular z axis preferentially in sample xy plane
Exp.Ordering = @(alpha,beta,gamma) f(beta,pi/2,w_beta);
[B,spc_xy] = pepper(Sys,Exp,Opt);

% (2) Partially oriented, molecular z axis preferentially along sample x axis
Exp.Ordering = @(alpha,beta,gamma) f(beta,pi/2,w_beta).*f(gamma,0,w_gamma);
[B,spc_x] = pepper(Sys,Exp,Opt);

% (3) Partially oriented, molecular z axis preferentially along sample y axis
Exp.Ordering = @(alpha,beta,gamma) f(beta,pi/2,w_beta).*f(gamma,pi/2,w_gamma);
[B,spc_y] = pepper(Sys,Exp,Opt);

% (4) Partially oriented, molecular z axis preferentially along sample z axis
Exp.Ordering = @(alpha,beta,gamma) f(beta,0,w_beta);
[B,spc_z] = pepper(Sys,Exp,Opt);

% Plotting
plot(B,spec_iso,B,spc_x,B,spc_y,B,spc_z,B,spc_xy);
xlabel('magnetic field (mT)')
legend('isotropic','ordered along x','ordered along y','ordered along z','ordered in xy plane');
legend boxoff
title('Partial ordering')
