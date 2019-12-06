% HYSCORE of a 14N nucleus (saffron)
%==========================================================================

clear, clc, clf

% Spin system with S=1/2 and one 14N nucleus
Sys.Nucs = '14N';
Sys.A_ = [0.8 0.5];  % aiso and T, MHz
Sys.Q = 0.4;         % e^2Qq/h, MHz

% HYSCORE experiment parameters
Exp.Field = 330;     % mT
Exp.Sequence = 'HYSCORE';
Exp.tau = 0.080;     % us
Exp.dt = 0.120;      % us
Exp.nPoints = 256;

% Increase number of orientations to get smooth powder average
Opt.nKnots = 300;    % 300 orientations over 90 degrees = 0.3 degrees

saffron(Sys,Exp,Opt);
