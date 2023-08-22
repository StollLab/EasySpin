% Single-crystal spectra of inter-system crossing spin-polarized triplet
%==========================================================================
clear, clf, clc
  
% Experimental parameters and simulation options
%------------------------------------------------------------------
Exp.mwFreq = 9.7;  % GHz
Exp.Range = [0 500];  % mT
Exp.Harmonic = 0;

% Orientation of the spin center in the crystal frame
Exp.MolFrame = [0 0 0];

Opt.Verbosity = 1;

% Defining the triplet spin system
%-------------------------------------------------------------------------------
Sys.S = 1;
Sys.g = 2;
Sys.lw = 5;  % mT
Sys.initState = {[0 1 0],'xyz'};  % selective population of the Ty state

D_cm = 0.06;  % cm^-1
D = unitconvert(D_cm,'cm^-1->MHz');
EoD = -0.2;  % E/D
Sys.D = D*[1 EoD];

% Simulations
%-------------------------------------------------------------------------------
% Orientation of the crystal frame in the lab frame.
% Change this value to turn the crystal in the spectrometer.
ang_z = deg2rad([0 0 0]);
ang_x = deg2rad([0 90 0]);
ang_y = deg2rad([0 90 90]);

Exp.SampleFrame = ang_z;
[B,spc_z] = pepper(Sys,Exp,Opt);
Exp.SampleFrame = ang_x;
[B,spc_x] = pepper(Sys,Exp,Opt);
Exp.SampleFrame = ang_y;
[B,spc_y] = pepper(Sys,Exp,Opt);

% Powder
Exp.MolFrame = [];
Exp.SampleFrame = [];
[B,spc_powder] = pepper(Sys,Exp,Opt);


% Prep for levelsplot
%-------------------------------------------------------------------------------
% phi/theta/chi angles for levelsplot are the Euler angles for the
% transformation from molecular lab frame to lab frame and are the inverse of
% Exp.SampleFrame (i.e. invert order and sign).
ptc_x = -fliplr(ang_x);  % only works for Exp.MolFrame = [0 0 0]
ptc_y = -fliplr(ang_y);  % only works for Exp.MolFrame = [0 0 0]
ptc_z = -fliplr(ang_z);  % only works for Exp.MolFrame = [0 0 0]

% Plotting
%-------------------------------------------------------------------------------
spc_crystal = {spc_x,spc_y,spc_z};
ptc = {ptc_x,ptc_y,ptc_z};
title_text = {'x','y','z'};
for k = 1:3
  subplot(2,3,k);
  levelsplot(Sys,ptc{k},Exp.Range,Exp.mwFreq);
  title(title_text{k});
  subplot(2,3,k+3);
  h = plot(B,spc_powder/max(spc_powder)*max(spc_crystal{k}),B,spc_crystal{k});
  h(2).Color = h(1).Color;
  h(1).Color = [1 1 1]*0.7;
  xlabel('magnetic field (mT)')
end
