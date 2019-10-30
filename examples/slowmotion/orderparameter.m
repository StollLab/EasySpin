% Calculate order parameter S(2,0)
%===============================================================================

% This example shows calculates an order parameter S(2,0) given an orientation
% potential coefficient (lambda_200 or C_200, i.e. L = 2, M = 0, and K = 0).

clear, clc, clf

% Set up orienational potential and population functions
%-------------------------------------------------------------------------------

% orientation potential coefficient and basis function
lambda20 = 2; % coefficient, unitless (multiples of kT)
Y20 = @(theta) (3*cos(theta).^2-1)*sqrt(5/pi)/4; % the associated orientational function

% set up the potential and the resulting probability distribution
T = 293; % temperature, K
kT = boltzm*293; % J
U = @(theta) -kT*lambda20*Y20(theta); % potential energy function
Z = integral(@(theta)exp(-U(theta)/kT).*sin(theta),0,pi/2); % partition function
P = @(theta) exp(-U(theta)/kT)/Z; % probability distribution (normalized)

% Calculation of ordering parameter for L=2, K=0
%-------------------------------------------------------------------------------
% To get order parameter S_(2,0), compute the integral
%                   \int_0^pi P(theta) Y20(theta) sin(theta) dtheta
% S_(2,0) = <Y20> = -----------------------------------------------
%                          \int_0^pi P(theta) sin(theta) dtheta
%
% The denominator integral is Z, and we already included it in P.

S20 = integral(@(theta)P(theta).*Y20(theta).*sin(theta),0,pi/2);


% Plotting
%-------------------------------------------------------------------------------
theta = linspace(0,pi/2,1001); % radians

subplot(2,1,1)
plot(theta*180/pi,P(theta));
xlabel('\theta (deg)');
ylabel('probability density');
ylim([0 max(P(theta))]);
title(sprintf('\\lambda_{20} = %g    S_{20} = %g',lambda20,S20));

subplot(2,1,2)
plot(theta*180/pi,U(theta)/echarge);
xlabel('\theta (deg)');
ylabel('potential energy (eV)');
