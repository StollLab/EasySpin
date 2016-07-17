% Susceptibility for S=1/2 using curry and explicit Curie law
%================================================================
% This example compares the results of a calculation with curry
% with the explicit form of the Curie law (which is accurate only
% in the high temperature region).

clear

% Quantities needed
S = 1/2;
g = 2;
B = 5000;  % mT
T = 2:0.5:100; % K

% Molar magnetic susceptibility according to Curie law
chi_Curie = mu0*avogadro*bmagn^2*g^2*S*(S+1)/3/boltzm./T;  % SI units

% Define the spin system and experimental parameters for curry
Sys.S = S;
Sys.g = g;
Exp.Field = B;
Exp.Temperature = T;

% Calculate magnetic susceptibility using curry
[muz,chizz] = curry(Sys,Exp);

% Plot the result. You can see that at low temperatures the Curie
% law is a bad approximation to the correct susceptibility.
plot(T,chizz,T,chi_Curie);
legend('exact','Curie law');
legend boxoff
xlabel('temperature (K)');
ylabel('molar magnetic susceptibility (m^3 mol^{-1})');
