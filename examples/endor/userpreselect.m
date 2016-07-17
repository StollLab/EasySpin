% User-supplied orientation pre-selection
%======================================================================
% This is an advanced example on how to speed up ENDOR simulations
% of multinuclear systems without using perturbation theory.
% It is based on a manual pre-computation of the orientation
% selectivity of the ENDOR experiment.

clear, clf

% Definition of the spin system
System.S = 1/2;
System.g = [2, 2, 2.2];
System.Nucs = '63Cu,1H';
System.A = [[40 40 400];  4+[-1 0 2]*1.5];
System.lwEndor = 0.02;

% Definition of the experiment
Experiment.mwFreq = 9.5;
Experiment.Range = [-1 1]*5 + 13.3;
Experiment.ExciteWidth = 50;
Experiment.Field = 310;

% ENDOR simulation options
Options.Verbosity = 1;
Options.nKnots = 91;
Options.Nuclei = 2;   % compute only 1H ENDOR

% ======= The important line ==============
Options.OriWeights = orisel(System,Experiment,Options);
% orisel(...) computes the orientation selectivity
% for the field/frequency values. The resulting weights
% are stored in Options.OriWeights, which are used by
% salt(...). This saves salt() the trouble of doing the
% orientation selection.

% In the following we only change the 1H coupling,
% which does not affect the orientation selection.
% As a consequence, salt() is much faster with
% Options.OriWeights given.

rho = linspace(0.01,1/3,10);
for k = 1:numel(rho)
  % change rhombicity of hyperfine coupling
  System.A(2,:) = 4 + [-1-rho(k) -1+rho(k) 2]*1.5;
  [x,y(k,:)] = salt(System,Experiment,Options);
end
y = y/max(y(:));  % normalize to amplitude 1

% Plotting
stackplot(x,y);

xlabel('frequency [MHz]');
title('ENDOR spectra');
