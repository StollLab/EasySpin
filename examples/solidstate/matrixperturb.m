% matrix diagonalization vs perturbation theory
%---------------------------------------------------------------

% In this example, the simulation of a powder EPR spectrum of
% a Cu2+ system with two nitrogen nuclei is done once with
% matrix diagonalization and then with second-order perturbation
% theory. The first method is more accurate, but the second
% is orders of magnitude faster.

clear, clf, clc

% Spin parameters
Sys.g = [2 2.2];
Sys.Nucs = '63Cu,14N,14N';
Sys.A = [50 300; 30 40; 30 40];
Sys.lwpp = 0.8;

% Experimental parameter
Exp.mwFreq = 10;
Exp.Range = [290 390];

% Two simulations with two methods (tic and toc measure the time)
Opt.Method = 'matrix';
Opt.Threshold = 0;
tic
[B,spc_matrix] = pepper(Sys,Exp,Opt);
toc

Opt.Method = 'perturb';
tic
[B,spc_perturb] = pepper(Sys,Exp,Opt);
toc

% Plotting
plot(B,spc_perturb,B,spc_matrix);
legend('perturb','matrix');
xlabel('magnetic field (mT)');
