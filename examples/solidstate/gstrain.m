% effect of g strain at various mw frequencies
%==========================================================================
clear, clf

% Spin system, experiment parameters and options
%------------------------------------------------------------
gStrain = [0.001 0.0008 0.0005];
Sys.g = [2.0104 2.0074 2.0026];
Exp.Harmonic = 0;

% Frequencies [GHz] and associated magnetic field ranges [mT]
%------------------------------------------------------------
Freqs = [3 9.5 35 95 350];
Ranges = [102 112; 334 344; 1240 1255; 3372 3395; 12425 12505];
nFreqs = numel(Freqs);

% Simulating all spectra with and without g strain
%------------------------------------------------------------
for k=1:nFreqs
  Exp.mwFreq = Freqs(k);
  Exp.Range = Ranges(k,:);
  
  Sys.gStrain = [0 0 0];
  [x{k},y1{k}] = pepper(Sys,Exp);
  Sys.gStrain = gStrain;
  [x{k},y2{k}] = pepper(Sys,Exp);
end

% Graphical rendering of the results
%------------------------------------------------------------
for k = 1:nFreqs
  subplot(nFreqs,1,k);
  h = plot(x{k},y1{k}/max(y1{k}),'r',x{k},y2{k}/max(y2{k}),'k');
  axis tight
  xx = xlim;
  text(xx(1),mean(ylim),sprintf('  %g GHz',Freqs(k)),'FontSize',8);
end
xlabel('magnetic field [mT]');
