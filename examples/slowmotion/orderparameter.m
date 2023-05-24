% Calculate order parameter S(2,0)
%===============================================================================

% This example calculates the order parameter S(2,0) from the orientational
% potential coefficient for (L,M,K) = (2,0,0) (as used in Sys.Potential
% when simulating slow-motion spectra using chili()).

% To get the order parameter S_(2,0), we need to compute the integral
%
%     S_(2,0) = <Y20> = \int_0^pi P(theta) Y20(theta) sin(theta) dtheta
%
% where Y20 is the orientational basis function with L=2, M=0, K=0.
% P(theta) is the orientational distribution given by
%
%     P(theta) = exp(-U(theta)/kT)/Z
%
% where Z is the integral
%
%     Z = \int_0^pi P(theta) sin(theta) dtheta

clear, clc, clf

% Set up orienational potential and population functions
%-------------------------------------------------------------------------------

% orientational potential coefficient for L=2, M=K=0
lambda200 = -200.1;  % in units of kB*T

% associated orientational basis function
Y20 = @(theta) (3*cos(theta).^2-1)/2;  % spherical harmonic (2,0)
D200 = Y20;  % Wigner D-function (2,0,0) is equal to spherical harmonic (2,0)

% orientational potential energy function
T = 293;  % temperature, K
kT = boltzm*T;  % J
U = @(theta) -kT*lambda200*D200(theta);

% partition function
Z = integral(@(theta)exp(-U(theta)/kT).*sin(theta),0,pi);

% normalized orientational distribution function
P = @(theta) exp(-U(theta)/kT)/Z;

% this integral should be equal to 1
P_integral = integral(@(theta)P(theta).*sin(theta),0,pi)

% Calculate the order parameter
%-------------------------------------------------------------------------------
f = @(theta) P(theta).*D200(theta).*sin(theta);
S20 = integral(f,0,pi);

% Plotting
%-------------------------------------------------------------------------------
theta = linspace(0,pi,181); % radians

subplot(2,1,1)
plot(theta*180/pi,P(theta));
xlabel('\theta (deg)');
ylabel('probability density');
ylim([0 max(P(theta))]);
title(sprintf('\\lambda_{2,0,0} = %g kT      S_{2,0,0} = %g',lambda200,S20));
grid on

subplot(2,1,2)
plot(theta*180/pi,U(theta)/echarge);
xlabel('\theta (deg)');
ylabel('potential energy (eV)');
grid on
