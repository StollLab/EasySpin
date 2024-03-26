% EPR spectrum of the S2 state of the oxygen evolving complex
% (Mn4CaO5) from photosystem II
%-------------------------------------------------------------
% This example demonstrates the recommended use of perturbation
% theory for the simulation of spectra from systems with many
% nuclei.

%55Mn hyperfine coupling tensors from
% Peloquin et al, J.Am.Chem.Soc. 2000, 122, 10926-10942
% https://doi.org/10.1021/ja002104f

clear

% 55Mn hyperfine couplings, all in MHz
AMn1 = -245 + [-1 -1 2]*(-13);
AMn2 = +217 + [-1 -1 2]*17;
AMn3 = -297 + [-1 -1 2]*14;
AMn4 = +200 + [-1 -1 2]*20;

Sys.g = [1.98 1.98 2.0];
Sys.lwpp = 2;
Sys = nucspinadd(Sys,'55Mn',AMn1);
Sys = nucspinadd(Sys,'55Mn',AMn2);
Sys = nucspinadd(Sys,'55Mn',AMn3);
Sys = nucspinadd(Sys,'55Mn',AMn4);

Exp.CenterSweep = [340 250];  % mT
Exp.mwFreq = 9.5;   % GHz

SimOpt.GridSize = 25;
SimOpt.Method ='perturb';

pepper(Sys,Exp,SimOpt);
