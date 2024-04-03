% A partially oriented Cu(II) complex
%==========================================================================
% This example shows how to use Exp.Ordering to specify partial ordering
% for a spin center.

clear, clf

% Paramagnetic complex
Sys.g = [2 2 2.2];
Sys.lwpp = 0.8;  % mT
Sys.Nucs = '63Cu';
Sys.A = [50 50 400];  % MHz

% Experimental parameters
Exp.mwFreq = 9.5;  % GHz
Exp.Range = [280 350];  % mT

% (1) Isotropic orientational distribution
[B,spc_iso] = pepper(Sys,Exp);

% (2) Partially oriented, molecular z axis predominantly in sample xy plane
Exp.Ordering = -2;
[B,spc_xy] = pepper(Sys,Exp);

% (3) Strongly oriented, molecular z axis predominantly along sample z axis
Exp.Ordering = 4;
[B,spc_z] = pepper(Sys,Exp);

% Plotting
plot(B,spc_iso,B,spc_xy,B,spc_z);
legend('isotropic distribution','preferential in xy plane', ...
  'preferential along z axis','location','best');
legend boxoff
xlabel('magnetic field (mT)')
