% Temperature dependence of magnetic moment and susceptibility
%================================================================
% This example shows how to use curry for the calculation
% of magnetic moment and magnetic susceptibility curves.

clear

% Define the spin system
Sys.S = 1/2;
Sys.g = 2;

% Define the experimental parameters: field and temperature
Exp.Field = 10000;
Exp.Temperature = 1:300;

% Call curry without outputs: curry will plot the results
curry(Sys,Exp);
