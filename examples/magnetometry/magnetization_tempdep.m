% Temperature-dependence of magnetization
%===============================================================

clear, clc, clf

% Very simple toy system with isotropic g factors
% and isotropic exchange
Sys.S = [1/2 1/2];
Sys.g = [2 2];
Sys.J = 5*30e3;

% A few temperature, and a field range
Exp.Temperature = [1 2 4 10];
Exp.Field = linspace(0,10,80)*1e3;

% calculate magnetization
mu_z = curry(Sys,Exp);

% plot
plot(Exp.Field/1e3,mu_z);
xlabel('magnetic field (T)');
ylabel('longitudinal magnetic moment (\mu_B)');
legend('1 K','2 K','4 K','10 K');
