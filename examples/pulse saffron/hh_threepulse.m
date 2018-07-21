% Three-pulse ESEEM of two 1H nuclei (saffron)
%==========================================================================

clear

Sys.Nucs = '1H,1H';
Sys.A = [2, -1; 5, 3]; % MHz

Exp.Sequence = '3pESEEM';
Exp.Field = 352.3; % mT
Exp.dt = 0.010; % increment in us
Exp.tau = 0.1; % delay between 1st and 2nd pulse in us
Exp.T = 0.06; % initial delay between 2nd and 3rd pulse, us

saffron(Sys,Exp);
