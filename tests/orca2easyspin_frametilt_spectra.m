function ok = test(opt)
% Check agreement of spectra simulated for ORCA calculations performed on
% identical structures but in two different ORCA coordinate frames
%----------------------------------------------------------------------------------------

% Methyl radical (tests g- and A-Frames)
%----------------------------------------------------------------------------------------

% Load file that contains a planar methyl radical
Sys = orca2easyspin('orca\methyl_all.oof');

% Load file containing the same radical but with C3 axis tilted away
% from the z axis of the ORCA frame.
Systilted = orca2easyspin('orca\methyl_tilted_all.oof');

Sys.lwpp = 0.1; % mT
Systilted.lwpp = Sys.lwpp;

Exp.mwFreq = 9.7; % GHz
Exp.Range = [335 360];
Exp.nPoints = 2e3;
Exp.Harmonic = 0;

[B,sim] = pepper(Sys,Exp);
[B,simtilted] = pepper(Systilted,Exp);

% Check that the EPR spectra match
ok(1) = areequal(sim,simtilted,0.1,'rel');

if opt.Display
  nexttile
  plot(B,sim,B,simtilted,B,sim-simtilted)
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  legend('sim','simtilted','residuals')
end

% Quinoline triplet (tests D- and A-Frames)
%----------------------------------------------------------------------------------------

% Load file that contains a quinoline triplet state
Sys = orca2easyspin('orca\quinoline_triplet.out');

% Load file containing the same triplet tilted and rotated
% in the ORCA frame.
Systilted = orca2easyspin('orca\quinoline_triplet_tilted.out');

Sys.lwpp = 0.1; % mT
Systilted.lwpp = Sys.lwpp;

Exp.mwFreq = 2; % GHz
Exp.Range = [0 150];
Exp.nPoints = 5e3;
Exp.Harmonic = 0;

[B,sim] = pepper(Sys,Exp);
[B,simtilted] = pepper(Systilted,Exp);

% Check that the EPR spectra match
ok(2) = areequal(sim,simtilted,0.1,'rel');

if opt.Display
  nexttile
  plot(B,sim,B,simtilted,B,sim-simtilted)
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  legend('sim','simtilted','residuals')
end