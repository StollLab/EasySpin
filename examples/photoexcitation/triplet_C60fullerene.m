% Exited triplet state EPR of C60 fullerene
%===========================================================
% see
%  Wasielewski et al, J.Am.Chem.Soc. 113, 2774-2776 (1991)
%     http://dx.doi.org/10.1021/ja00007a074
%  Dauw et al, J.Chem.Phys. 112, 7102-7110 (2000)
%     http://dx.doi.org/10.1063/1.481305

clear

Sys.S = 1;
Sys.g = 2.00; % neglecting small anisotropy

% Wasielewski 1991 values for zero-field splitting
D = 0.0114 * 100*clight/1e6;     % cm^-1 -> MHz
E = 0.00069 * 100*clight/1e6;    % cm^-1 -> MHz

Sys.D = [D E]; % MHz

Sys.lwpp = 0.9;   % mT

Exp.mwFreq = 9.248;  % GHz
Exp.Range = [310 350];      % mT

% The three numbers in Exp.Temperature give the relative
% populations of the three sublevels at zero-field (where
% T0 is lowest in energy) in order of increasing energy.
Exp.Temperature = [1 0 0];

pepper(Sys,Exp);
