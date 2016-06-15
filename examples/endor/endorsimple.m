% 1H ENDOR of a powder
%=====================================================================
clear, clf

% Spin system
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A_ = [1 3];
Sys.lwEndor = 0.3;

% Experimental parameters
Exp.Field = 350.1;
% frequency range is determined automatically

% Simulation
salt(Sys,Exp);
