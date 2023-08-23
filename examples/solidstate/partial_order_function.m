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

% (1) Isotropic orientational distribution
[B,spec_iso] = pepper(Sys,Exp,Opt);

% Distribution widths along orientational angles
% (beta and gamma are the second and third of the Euler angles that
% transform the sample frame into the molecular frame)
w_beta = deg2rad(30);
w_gamma = deg2rad(30);

% Elementary distribution function for an angle x
gauss = @(x,w)exp(-(x./w).^2);

% (1) Partially oriented, B0 preferentially in molecular xy plane
beta0 = -pi/2;
Exp.Ordering = @(beta) gauss(beta-beta0,w_beta);
[B,spec_xy] = pepper(Sys,Exp,Opt);

% (2) Partially oriented, B0 preferentially along molecular x axis
Exp.Ordering = @(beta,gamma) gauss(beta-beta0,w_beta).*gauss(gamma,w_gamma);
[B,spec_x] = pepper(Sys,Exp,Opt);

% (3) Partially oriented, B0 preferentially along molecular y axis
gamma0 = -pi/2;
Exp.Ordering = @(beta,gamma) gauss(beta-beta0,w_beta).*gauss(gamma-gamma0,w_gamma);
[B,spec_y] = pepper(Sys,Exp,Opt);

% (4) Partially oriented, B0 preferentially along molecular z axis
Exp.Ordering = @(beta) gauss(beta,w_beta);
[B,spec_z] = pepper(Sys,Exp,Opt);

% Plotting
plot(B,spec_iso,B,spec_x,B,spec_y,B,spec_z,B,spec_xy);
legend boxoff
xlabel('magnetic field (mT)')
legend('isotropic','along x','along y','along z','in xy plane');
