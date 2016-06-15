% matrix diagonalization vs perturbation theory
%---------------------------------------------------------------

% In this example, the simulation of a powder EPR spectrum of
% a Cu2+ system with two nitrogen nuclei is done once with
% matrix diagonalization and then with second-order perturbation
% theory. The first method is more accurate, but the second
% is orders of magnitude faster.

clear, clf

% spin parameters
Sys.g = [2 2.2];
Sys.Nucs = '63Cu,14N,14N';
Sys.A = [50 300; 30 40; 30 40];
Sys.lwpp = 0.8;

% experimental parameter
Exp.mwFreq = 10;
Exp.Range = [290 390];

% two simulations with two methods (tic and toc measure the time)
Opt.Method = 'matrix';
Opt.Threshold = 0;
tic
[x,y1] = pepper(Sys,Exp,Opt);
toc

Opt.Method = 'perturb';
tic
[x,y2] = pepper(Sys,Exp,Opt);
toc

% graphical plotting
plot(x,y2,'r',x,y1,'b');
