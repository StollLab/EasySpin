% Slow-motional EPR spectrum of TEMPONE
%===============================================================================
% This reproduces Figure 5 on page 38 of
%   D.J.Schneider, J.H.Freed
%   Calculating Slow Motional Magnetic Resonance Spectra
%   Biol.Magn.Reson. 8 (1989), 1-76
% (for parameters, see Table 2 on page 39, No. 2)

clear, clf

Tempone.g = [2.0088 2.0061 2.0027];
Tempone.Nucs = '14N';
Tempone.A = mt2mhz([5.8,5.8,30.8]/10); % G -> mT -> MHz

Tempone.Diff = 1e6; % s^-1
Tempone.lw = [0 0.2]; % Lorentzian broadening, FWHM, mT

Experiment.mwFreq = 9.2646; % GHz

Opt.LLKM = [14 7 6 2];

[B,spc] = chili(Tempone,Experiment,Opt);

% Plotting
%-------------------------------------------------------------------------------
B0 = 330; % mT
plot(B-B0,spc);
axis tight
xlabel('magnetic field (mT)');
