% Slow-motional EPR spectrum of TEMPONE
%=========================================================
% This reproduces Figure 5 on page 38 of
% D.J.Schneider, J.H.Freed
% Calculating Slow Motional Magnetic Resonance Spectra
% Biol.Magn.Reson. 8 (1989), 1-76

clear, clf

Tempone.g = [2.0088 2.0061 2.0027];
Tempone.Nucs = '14N';
Tempone.A = mt2mhz([5.8,5.8,30.8]/10);

Tempone.Diff = 1e6;
Tempone.lw = 0.2;

Experiment.mwFreq = 9.2646;

chili(Tempone,Experiment);
