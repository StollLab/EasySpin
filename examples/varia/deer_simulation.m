% Simulate a 4-pulse DEER signal
%===============================================================================
clear, clc

% Set up time and distance vectors
t = (-0.1:0.01:3).';  % dipolar evolution time, µs
nr = numel(t);
r = linspace(1.5,8,nr);  % spin-spin distance, nm

% Build distance distribution
r0 = 3.5;  % center of Gaussian, nm
fwhm = 0.6;  % full width at half maximum, nm
P = gaussian(r,r0,fwhm).';

% Calculate oscillating time trace
K = dipkernel(t,r);  % dipolar kernel; K(:,i) contains time signal for distance r(i)
V_intra0 = K*P;

% Add modulation depth and background
lam = 0.3;  % modulation depth
tauB = 6;  % decay time of the expinential background decay, µs
V_inter = exp(-abs(t)/tauB);
V = ((1-lam)+lam*V_intra0).*V_inter;

% Plotting
subplot(2,1,1)
plot(t,V,t,(1-lam)*V_inter,'--');
xlabel('time (µs)')
ylabel('echo amplitude (norm.)')
grid on
title('DEER time trace')

subplot(2,1,2)
plot(r,P)
xlabel('distance (nm)')
ylabel('density (nm^{-1})')
grid on
title('Distance distribution')
