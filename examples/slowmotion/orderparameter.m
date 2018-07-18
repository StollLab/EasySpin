% Calculate order parameter S(2,0)
%===============================================================================

% How to calculate an order parameter S(2,0) given an orientation potential
% coefficient (lambda_20 or C_20, i.e. L = 2 and K = 0; M = 0).

clear, clc, clf

% Set up orienational potential and population functions
%-------------------------------------------------------------------------------
theta = linspace(0,pi/2,1001); % radians
ct = cos(theta);

% orientations (pseudo)potential coefficient
C20 = 5; % unitless
Y20 = (3*ct.^2-1)*sqrt(5/pi)/4; % the associated orientational function

% set up the potential and the resulting probability distribution
T = 293; % temperature, K
kT = boltzm*293; % J
U = -kT*C20*Y20; % potential energy function
P = exp(-U/kT); % probability distribution (unnormalized)

% normalize the probability distribution
dtheta = theta(2)-theta(1);
Z = trapz(P.*sin(theta))*dtheta; % partition function
P = P/Z;

% Calculation of ordering parameter for L=2, K=0
%-------------------------------------------------------------------------------
% basis function for axial order: Y_(2,0)
% (is the simplest possible; can be normalized or not)
Y20 = (3*ct.^2-1)*sqrt(5/pi)/4;

% To get order parameter S_(2,0), compute the integral
%                   \int_0^pi P(theta) Y20(theta) sin(theta) dtheta
% S_(2,0) = <Y20> = -----------------------------------------------
%                          \int_0^pi P(theta) sin(theta) dtheta

S20 = trapz(P.*Y20.*sin(theta))*dtheta;

% Plotting
%-------------------------------------------------------------------------------
subplot(2,1,1)
plot(theta*180/pi,U/echarge);
set(gca,'FontSize',12);
xlabel('theta (deg)');
ylabel('potential energy (eV)');

subplot(2,1,2)
plot(theta*180/pi,P);
set(gca,'FontSize',12);
xlabel('theta (deg)');
ylabel('probability density');
ylim([0 max(P)]);
legend(sprintf('S20 = %g',S20));
