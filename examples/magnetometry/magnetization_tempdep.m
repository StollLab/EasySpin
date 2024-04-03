% Temperature-dependence of magnetization
%===============================================================

clear, clc, clf

% Simple two-spin system with isotropic g factors and isotropic exchange
Sys.S = [1/2, 1/2];
Sys.g = [2, 2];
Sys.J = 5*30e3;  % cm^-1 -> MHz

% Define few temperature, and a field range
Exp.Temperature = [1 2 4 10];  % K
B = linspace(0,10,80);  % T
Exp.Field = B*1e3;  % mT

% Calculate magnetization
mu_z = curry(Sys,Exp);

% Plot
plot(B,mu_z);
xlabel('magnetic field (T)');
ylabel('longitudinal magnetic moment (\mu_B)');
legend('1 K','2 K','4 K','10 K');
grid on
