% Orientation-selective ENDOR of many protons
%=======================================================
% When many ENDOR nuclei are present in a spin system,
% then it is advisable to use perturbation theory, since
% matrix diagonalization is very slow.
clear, clf

% Spin system: several protons with random hyperfine tensorss
Sys.g = [2 2.1 2.2];
Sys.Nucs = '1H,1H,1H,1H,1H';
Sys.A = rand(5,3)*16;
Sys.AFrame = rand(5,3)*pi;
Sys.lwEndor = 0.1;

% Experiment
Exp.Field = 3280;
Exp.CenterSweep = [140 20];
Exp.mwFreq = 94;
Exp.ExciteWidth = 200;

% Simulation options
% - use perturbation theory instead of matrix diagonalization
Opt.Method = 'perturb1';
% - keep transitions separate
Opt.Output = 'separate';
% - high orientational resolution, no interpolation
Opt.nKnots = [61 0];

salt(Sys,Exp,Opt);
