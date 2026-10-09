% Mn(II) at W band with D strain: matrix diagonalization vs perturbation theory
%==========================================================================
% At high field, the zero-field splitting of Mn(II) is small compared to the
% electron Zeeman interaction, and second-order perturbation theory gives
% accurate spectra much faster than matrix diagonalization. The correlated
% D/E strain broadens the outer fine-structure transitions, whereas the
% central mS = -1/2 <-> +1/2 hyperfine sextet stays narrow, since it is
% affected by D only in second order. Perturbation theory computes only
% allowed transitions, so it misses the weak forbidden-transition doublets
% between the sextet lines.

clear, clf, clc

% Spin system
Sys.S = 5/2;
Sys.g = 2.0;
Sys.Nucs = '55Mn';
Sys.A = -250;  % MHz
Sys.D = [300 60];  % D and E, MHz
Sys.StrainPars = {'D(1)','D(2)'};  % strains of D and E
Sys.StrainFWHM = [150 50];  % MHz
Sys.StrainCorr = 0.3;  % correlation between D and E strains
Sys.lwpp = 0.5;  % mT

% Experimental parameters
Exp.mwFreq = 95;  % GHz
Exp.Range = [3250 3540];  % mT

% Simulations with two methods (tic and toc measure the time)
Opt.Method = 'matrix';
tic
[B,spc_matrix] = pepper(Sys,Exp,Opt);
toc

Opt.Method = 'perturb';
tic
[B,spc_perturb] = pepper(Sys,Exp,Opt);
toc

% Plotting
% The lower panel is zoomed vertically to show the broad outer transitions.
subplot(2,1,1);
plot(B,spc_matrix,B,spc_perturb);
axis tight

subplot(2,1,2);
plot(B,spc_matrix,B,spc_perturb);
axis tight
ylim([-1 1]*0.02*max(abs(spc_matrix)));
legend('matrix','perturb');
legend boxoff
xlabel('magnetic field (mT)');
