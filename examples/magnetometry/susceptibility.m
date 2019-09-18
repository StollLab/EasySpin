% Temperature dependence of magnetic susceptibility
%=========================================================

clear, clc, clf

Sys.S = 7/2;
D = 2; % cm^-1
Sys.D = D*30e3; % cm^-1 -> MHz

Exp.Field = 1000; % mT
T = [0.1:0.1:1, 1:0.5:50]; % K
Exp.Temperature = T;

% Calculate molar susceptibility
Opt.Output = 'chimol';
chimol = curry(Sys,Exp,Opt);

% Plot chi and chi*T
subplot(2,1,1)
plot(T,chimol);
xlabel('temperature (K)');
ylabel('\chi_{mol} (m^3 mol^{-1})');

subplot(2,1,2)
plot(T,chimol.*T);
xlabel('temperature (K)');
ylabel('\chi_{mol}T (K m^3 mol^{-1})');
