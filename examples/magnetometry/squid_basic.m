% Temperature dependence of magnetic moment and susceptibility
%================================================================
% This example shows how to use curry() for the calculation
% of magnetic moment and magnetic susceptibility curves.

clear

% Define the spin system
Sys.S = 1/2;
Sys.g = 2;

% Define the experimental parameters: field and temperature
Exp.Field = 10000;  % mT
Exp.Temperature = 1:300;  % K

% Call curry without requesting outputs: curry plots the results
curry(Sys,Exp);

%%
% repeat calculations using CGS units
Opt.Units = 'CGS';
curry(Sys,Exp,Opt);
