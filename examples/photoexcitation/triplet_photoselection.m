% Photoselection for spin-polarized triplet EPR spectrum
%===============================================================================

% This example demonstrates how to take into account different laser
% excitation geometries when simulating the transient EPR spectrum of a
% photo-generated spin-polarized triplet radical.

clear, clc, clf

% Define spin system
Triplet.S = 1;
Triplet.D = [900 -160];  % MHz
Triplet.lwpp = 1;  % mT
Triplet.tdm = 'x';  % orientation of tdm in molecular frame
Triplet.initState = {[0 0 1],'xyz'};  % non-equilibrium populations

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [280 400];  % mT
Exp.Harmonic = 0;  % for transient EPR, turn off field modulation

Opt.GridSymmetry = 'D2h';

% Simulate spectra under various forms of photoexcitation
Exp.lightBeam = '';
[B,spc_complete] = pepper(Triplet,Exp,Opt);

Exp.lightScatter = 0.2;  % include 20% isotropic contribution
Exp.lightBeam = 'perpendicular';
[B,spc_perp] = pepper(Triplet,Exp,Opt);
Exp.lightBeam = 'parallel';
[B,spc_para] = pepper(Triplet,Exp,Opt);
Exp.lightBeam = 'unpolarized';
[B,spc_unpol] = pepper(Triplet,Exp,Opt);

% Plotting
plot(B,spc_complete,B,spc_perp,B,spc_para,B,spc_unpol);
xlabel('magnetic field (mT)')
legend('complete','perpendicular','parallel','unpolarized')
legend boxoff
