% relaxation (spidyan)
%==========================================================================
% demonstration of how to switch relaxation on for specific events and how
% to set a userdefined equilibrium state
% relaxation types can be switched off, by removing them from Sys, or by
% setting to zero, e.g. to use only T1: Sys.T1 = 2; Sys.T2 = 0;

clear

% Default Spin System
Sys.S = 1/2;
Sys.ZeemanFreq = 9.500; % GHz
Sys.T1 = 2; % us
Sys.T2 = 1; % us
Sys.eqState = -sop(Sys.S,'z'); % make up some equilibrium state

% Pulse Definitions
Rectangular.Type = 'rectangular';
Rectangular.tp = 0.02;
Rectangular.Flip = pi;

% A default Experiment/Sequence
Exp.mwFreq = 9.5; % GHz
Exp.Sequence = {Rectangular 5};
Exp.DetOperator = {'z1'};

% Options
Opt.Relaxation = [0 1]; % switches relaxation on only during the free evolution period 

[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% plotting
figure(1)
clf
plot(TimeAxis,real(Signal))
xlabel('t (\mus)')
ylabel('<S_i>')
legend(Exp.DetOperator)