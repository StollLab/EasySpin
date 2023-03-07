% Single-crystal spectra of inter-system crossing spin-polarized triplet
%==========================================================================
clear, clf, clc
  
% Experimental parameters and simulation options
%------------------------------------------------------------------
Exp.mwFreq = 9.7;
Exp.Range = [0 450];
Exp.Harmonic = 0;

% Orientation of the crystal frame in the lab frame.
% Change this value to turn the crystal in the spectrometer.
Exp.SampleFrame = [0 0 0]*pi/180;

Opt.Verbosity = 0;

% Defining two spin systems, one with axial and one with rhombic D
%-----------------------------------------------------------------------
Sys1.S = 1;
Sys1.g = 2;
Sys1.lw = 5;  % mT
Sys1.initState = {[0 1 0],'xyz'};  % selective population ot the Ty state
Sys2 = Sys1;

D_cm = 0.06;  % cm^-1
D = D_cm*100*clight/1e6;  % cm^-1 -> MHz
Sys1.D = D*[1 0];
Sys2.D = D*[1 -0.1];

% phi/theta/chi angles for levelsplot are the Euler angles for the
% transformation from molecular lab frame to lab frame and are the inverse of
% Exp.SampleFrame.
ang = -fliplr(Exp.SampleFrame);

% Simulations & plotting
%-------------------------------------------------------------
subplot(3,2,1);
levelsplot(Sys1,ang,Exp.Range,Exp.mwFreq);
title('Axial D tensor');

subplot(3,2,3);
pepper(Sys1,Exp,Opt);

subplot(3,2,2);
levelsplot(Sys2,ang,Exp.Range,Exp.mwFreq);
title('Non-axial D tensor');

subplot(3,2,4);
pepper(Sys2,Exp,Opt);

Exp.SampleFrame = [];  % powder spectrum

subplot(3,2,5);
pepper(Sys1,Exp,Opt);

subplot(3,2,6);
pepper(Sys2,Exp,Opt);
