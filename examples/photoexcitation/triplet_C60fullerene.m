% Exited triplet state EPR of C60 fullerene
%===========================================================
% see
%  Wasielewski et al, J.Am.Chem.Soc. 113, 2774-2776 (1991)
%     https://doi.org/10.1021/ja00007a074
%  Dauw et al, J.Chem.Phys. 112, 7102-7110 (2000)
%     https://doi.org/10.1063/1.481305

clear

Sys.S = 1;
Sys.g = 2.00;  % neglecting small anisotropy

% Wasielewski 1991 values for zero-field splitting
D_cm = -0.0114;  % cm^-1
E_cm =  0.00069;  % cm^-1

Sys.D = unitconvert([D_cm E_cm],'cm^-1->MHz');  % cm^-1 -> MHz

Sys.lwpp = 0.9;  % mT

Exp.mwFreq = 9.248;  % GHz
Exp.Range = [310 350];  % mT

% The three numbers in the first input to Sys.initState give the
% relative populations of the three sublevels Tx, Ty and Tz
Sys.initState = {[0.5 0.5 0],'xyz'};

pepper(Sys,Exp);
