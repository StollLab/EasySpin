% Comparison between fast- and slow-motion simulation algorithms
%===============================================================================
clear, clf, clc

% Parameters
%-------------------------------------------------------------------------------
Nx.g = [2.0088, 2.0061, 2.0027];
Nx.Nucs = '14N';
Nx.A = unitconvert([5.8, 5.8, 30.8]/10,'mT->MHz');   % G -> MHz;

Exp.mwFreq = 9.5;    % GHz
Exp.CenterSweep = [338, 10];  % mT

tcorr = 1e-9;  % seconds
lw = 0.03;     % mT

% Simulation
%-------------------------------------------------------------------------------
% using chili (stochastic Liouville equation solver)
Nx.tcorr = tcorr;
Nx.lw = [0 lw]; % Lorentzian broadening
[B,spc_s] = chili(Nx,Exp);

% using garlic (perturbational treatment valid only in fast-motion limit)
Nx.lw = [0 lw];  % Lorentzian broadening
Nx.tcorr = tcorr;
[B,spc_f] = garlic(Nx,Exp);

% Plotting
%-------------------------------------------------------------------------------
plot(B,spc_s/max(spc_s),'k',B,spc_f/max(spc_f),'r');
xlabel('magnetic field [mT]');
legend('chili','garlic');
axis tight
title(sprintf('correlation time: %g ns',Nx.tcorr/1e-9));
