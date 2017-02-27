% [Cu(mnt)2] at 1.4 GHz, single-crystal
%==========================================================================
% Values from Kirmse et al., Inorg.Chem. 23, 3333-3338 (1984)
clear, clf

% spin system with Cu (mixture of 63Cu and 65Cu)
Sys.Nucs = 'Cu';
Sys.lwpp = 0.5;
Sys.g = [2.020 2.023 2.089];
A = [-43.1 -43.7 -171.4]*1e-4; % cm^-1
Sys.A = A*100*clight/1e6; % conversion cm^-1 -> MHz

% experimental parameters
Exp.mwFreq = 1.4;
Exp.Range = [10 90];

% options
Opt.Verbosity = 1;

% simulation of two single-crystal spectra
% in one pepper call
Exp.CrystalOrientation = [0 0 0; 0 pi/2 0];
Opt.Output = 'separate';
[B,spec] = pepper(Sys,Exp,Opt);
spec1 = spec(1,:);
spec2 = spec(2,:);

% display
plot(B,spec1,'r',B,spec2,'b'); axis tight
legend('0^o','90^o');
xlabel('magnetic field [mT]');
title(sprintf('%g GHz single-crystal spectra of Cu(mnt)2',Exp.mwFreq));
