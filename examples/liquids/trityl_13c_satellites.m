% Isotropic spectrum of trityl radical including 13C satellite lines
%==========================================================================
% Zoom in to see the 13C satellite lines!

% W. Moore et al,
% J. Magn. Reson. 318, 106797 (2020)
% https://doi.org/10.1016/j.jmr.2020.106797

clear, clc

Trityl.g = [2.0033 2.0032 2.0075];
Trityl.Nucs = 'C';
Trityl.A = [18 162];  % MHz
Trityl.lwpp = 0.1;  % mT

Exp.mwFreq = 9.6;  % GHz
Exp.Range = [340 345];  % mT

Opt.separate = 'components';

subplot(2,1,1)
garlic(Trityl,Exp,Opt);
title('fluid solution')

subplot(2,1,2)
pepper(Trityl,Exp,Opt);
title('frozen solution')
