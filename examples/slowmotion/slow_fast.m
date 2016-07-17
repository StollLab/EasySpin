% Comparison between fast- and slow-motion simulation algorithms
%==========================================================================
clear, clf

% Parameters
%----------------------------------------------------
g = [2.0088, 2.0061, 2.0027];
A = mt2mhz([5.8, 5.8, 30.8]/10);   % G -> MHz
Nx = struct('g',g,'Nucs','14N','A',A);

Par = struct('mwFreq',9.5,'CenterSweep',[338, 10]);

tcorr = 1e-9;  % seconds
lw = 0.05;      % mT

% Simulation
%----------------------------------------------------
% Slow-motion
Nx.tcorr = tcorr;
Nx.lw = lw; % mT -> MHz
[B,spc_s] = chili(Nx,Par);

% Fast-motion
Nx.lw = [0 lw];  % Lorentzian broadening
Nx.tcorr = tcorr;
[B,spc_f] = garlic(Nx,Par);

% Graphical rendering
%----------------------------------------------------
plot(B,spc_s/max(spc_s),'k',B,spc_f/max(spc_f),'r');
xlabel('magnetic field [mT]');
legend('chili','garlic');
axis tight
title(sprintf('correlation time: %g ns',Nx.tcorr/1e-9));
