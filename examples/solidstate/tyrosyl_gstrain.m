% Tyrosyl radical: effect of g strain at various mw frequencies
%==========================================================================
% The tyrosyl radical Y122 of E. coli ribonucleotide reductase has a small
% g anisotropy. A distribution of g values (g strain), largest for gx due
% to variations in hydrogen bonding, broadens the spectrum proportionally
% to the microwave frequency, so its effect becomes visible only at high
% frequencies. The hyperfine couplings are modeled as a frequency-independent
% Gaussian broadening (HStrain), which dominates at low frequencies.
%
% g values: Gerfen et al., J. Am. Chem. Soc. 115, 6420 (1993),
%   https://doi.org/10.1021/ja00067a071

clear, clf

% Spin system, experiment parameters and options
%------------------------------------------------------------
gFWHM = [0.001 0.0005 0.0003];  % FWHM of g distributions
Sys.g = [2.0091 2.0046 2.0022];
Sys.HStrain = [1 1 1]*40;  % MHz, unresolved hyperfine couplings
Exp.Harmonic = 0;

% Frequencies [GHz] and associated magnetic field ranges [mT]
%------------------------------------------------------------
Freqs = [3 9.5 35 95 263];
Ranges = [104 110; 335 342; 1242 1252; 3374 3394; 9345 9395];
nFreqs = numel(Freqs);

% Simulating all spectra with and without g strain
%------------------------------------------------------------
for k = 1:nFreqs
  Exp.mwFreq = Freqs(k);
  Exp.Range = Ranges(k,:);
  
  [B{k},spc1{k}] = pepper(Sys,Exp);
  SysStrain = Sys;
  SysStrain.StrainPars = {'g(1)','g(2)','g(3)'};
  SysStrain.StrainFWHM = gFWHM;
  [B{k},spc2{k}] = pepper(SysStrain,Exp);
end

% Graphical rendering of the results
%------------------------------------------------------------
for k = 1:nFreqs
  subplot(nFreqs,1,k);
  h = plot(B{k},spc1{k}/max(spc1{k}),'r',B{k},spc2{k}/max(spc2{k}),'k');
  axis tight
  xx = xlim;
  text(xx(1),0.8,sprintf('  %g GHz',Freqs(k)),'FontSize',8);
end
xlabel('magnetic field (mT)');
