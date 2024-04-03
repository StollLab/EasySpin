% Comparison of magnetic moment and effective magnetic moment
%===============================================================================
% This example illustrates the difference between the single-molecule magnetic
% moment and the effective magnetic moment. The first has SI units of J/T, and
% the latter is a unitless quantity.
%
% curry can calculate both. The magnetic moment can be calculated using
% Opt.Output = 'mu' (in J/T), or Opt.Output = 'muBM' (in units of Bohr
% magnetons). The effective magnetic moment is calculated by Opt.Output =
% 'mueff'.

clear, clc, clf

% Definition of the spin system
Sys.S = 1;
Sys.g = 2.05;
Sys.D = 5*30e4;  % cm^-1 -> MHz

% Field and temperature of the measurement
B = 0.1;  % T
T = 0:100;  % K
Exp.Field = B*1e3;  % mT
Exp.Temperature = T;

% Calculate both magnetic moment and effective magnetic moment with one
% call to curry.
Opt.Output = 'mu mueff';
[mu,mueff] = curry(Sys,Exp,Opt);

% The high-temperature limit of the effective magnetic moment is
% g*sqrt(S*(S+1)), assuming spin-only.
mueff_lim = Sys.g*sqrt(Sys.S*(Sys.S+1));

% Plotting
subplot(2,1,1);
plot(T,mu/bmagn);
xlabel('temperature (K)');
ylabel('\mu (\mu_B)');
title('magnetic moment');
grid on

subplot(2,1,2);
plot(T,mueff);
line([0 max(T)],mueff_lim*[1 1],'Color','k','LineStyle','--');
xlabel('temperature (K)');
ylabel('\mu_{eff}');
title('effective magnetic moment');
grid on

