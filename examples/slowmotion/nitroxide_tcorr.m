% Correlation time dependence of a nitroxide spectrum
%==========================================================================
clear, clf

% Parameters
%--------------------------------------------------------------------------
g = [2.008, 2.006, 2.003];
Nitroxide = struct('g', g, 'Nucs', '14N', 'A', [20, 20, 85]);
Nitroxide.lw = 0.2;

Experiment = struct('mwFreq',9.5);

Options = struct('Verbosity',1);
Options.LLKM = [24 14 6 6];

% Simulations
%--------------------------------------------------------------------------
tcorr = [0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000]*1e-9;
nSpc = numel(tcorr);
for k = 1:nSpc
  Nitroxide.tcorr = tcorr(k);
  [B{k},spc{k}] = chili(Nitroxide,Experiment,Options);
end

% Graphical rendering
%--------------------------------------------------------------------------
hold on

for k = 1:nSpc
  spcnorm = spc{k}/max(spc{k});
  shift = k-1;
  plot(B{k},spcnorm + shift);
end

axis tight
xlabel('magnetic field [mT]');
