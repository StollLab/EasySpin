% Temperature dependence of magnetic susceptibility
%=========================================================

clear, clc, clf

Sys.S = 7/2;
Sys.D = 2*30e3; % MHz

Exp.Field = 1000; % mT
T = [0.1:0.1:1, 1:0.5:50]; % K
Exp.Temperature = T;

% Calculate magnetic moment mu and susceptibility chi
[mu,chi] = curry(Sys,Exp);

% Plot chi, chi*T, and 1/chi
subplot(3,1,1)
plot(T,chi);
xlabel('temperature (K)');
ylabel('\chi_{mol} (m^3 mol^{-1})');

subplot(3,1,2)
plot(T,chi.*T);
xlabel('temperature (K)');
ylabel('\chi_{mol}T (K m^3 mol^{-1})');

subplot(3,1,3)
plot(T,1./chi);
xlabel('temperature (K)');
ylabel('1/\chi_{mol} (mol m^{-3})');
