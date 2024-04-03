% Temperature dependence of magnetic moment and susceptibility
%================================================================
% This example shows calculates magnetic moment and magnetic
% susceptibility curves for a system of two coupled spins-1/2.
%
% It shows how to deal with exchange coupling definition different
% from EasySpin, and how to convert units.

clear, clf, clc

% Define the spin system
Sys.S = [1/2, 1/2];
Sys.g = [2, 2];

% Here is a J value, given for -2*J12*S1*S2, and in units of cm^-1
J12 = 30;   % cm^-1, H = -2*J12*S1*S2

% First, we convert units to the required MHz
J12 = unitconvert(J12,'cm^-1->MHz');  % cm^-1 -> MHz

% Then, we set the exchange coupling. Since EasySpin defines the
% exchange operator as H = J*S1*S2, we set J = -2*J12.
Sys.J = -2*J12;

% Define the experimental parameters: field and temperature
Exp.Field = linspace(0,10,50)*1000;  % mT
Exp.Temperature = 4;   % K

% Call curry and plot the result
[muz,chizz] = curry(Sys,Exp);

% Plot the results
B = Exp.Field/1e3;  % mT -> T

subplot(1,2,1);
plot(B,muz);
xlabel('magnetic field (T)');
ylabel('magnetic moment (\mu_B)');
grid on

subplot(1,2,2);
plot(B,chizz);
xlabel('magnetic field (T)');
ylabel('magnetic susceptibility (m^3 mol^{-1})');
grid on
