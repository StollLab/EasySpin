% Calculate order parameter S(2,0)
%==========================================================================

% How to calculate an order parameter S(2,0)
% given an ordering potential coefficient
% (lambda_20 or C_20)

clear, clc, clf

theta = linspace(0,pi/2,1001);

% ordering (pseudo)potential coefficient
C20 = -50;
Y20 = (3*cos(theta).^2-1)/2;

% set up the potential and the resulting
% probability distribution
kT = boltzm*293;
U = -kT*C20*Y20;
P = exp(-U/kT);

% normalize the probability distribution
dtheta = theta(2)-theta(1);
P = P/sum(P.*sin(theta))/dtheta;

% plotting
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

%%
% basis function for axial order: Y_(2,0)
% (is the simplest possible; has to be
% normalized)
Y20 = (3*cos(theta).^2-1)/2;

% To get order parameter S_(2,0), compute
% integral S_(2,0) = int_0^pi P(theta)*Y20(theta) sin(theta) dtheta

S20 = sum(P.*Y20.*sin(theta))*dtheta
