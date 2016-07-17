% Bridged biaryl cation radical (1000s of lines)
%==========================================================================
% Cation radical of 6-hydrodipyrido[1,2-c:2',1'-e]-imidazole
% The spectrum contains a total of 9216 resonance lines.

clear, clf

% Hyperfine couplings in units of Gauss (from literature)
A_Gauss = [+4.34,-2.39,-0.65,-2.81,-0.23,+24.24];
A_MHz = mt2mhz(A_Gauss/10,gfree); % xonversion to MHz

% The spin system contains a total of 12 nuclei, in groups of 2
Sys.g = 2.00316;
Sys.Nucs = '14N,1H,1H,1H,1H,1H';
Sys.n = [2 2 2 2 2 2];
Sys.A = A_MHz;
Sys.lwpp = [0,0.006];

Exp.mwFreq = 9.532;
Exp.nPoints = 8192;

% Simulation and plotting of the solution cw EPR spectrum
garlic(Sys,Exp);
