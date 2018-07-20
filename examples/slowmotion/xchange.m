% Effect of Heisenberg exchange in a nitroxide sample
%===============================================================================
clear, clf, clc

% Nitroxide and experimental parameters
%-------------------------------------------------------------------------------
Nx.g = [2.0088, 2.0061, 2.0027];
Nx.Nucs = '14N';
Nx.A = [16 16 86];   % MHz
Nx.tcorr = 0.2e-9;   % ns
Nx.lw = 0.3;         % mT

Exp.mwFreq = 9.5;  % GHz
Exp.CenterSweep = [338, 10];  % mT

% Loop over range of exchange frequencies, simulation
%-------------------------------------------------------------------------------
ExchangeFreq = 10.^linspace(1,3,10); % in MHz
for k = 1:numel(ExchangeFreq)
  Nx.Exchange = ExchangeFreq(k);
  [B,spc(k,:)] = chili(Nx,Exp);
end

% Graphical rendering
%-------------------------------------------------------------------------------
stackplot(B,spc);
xlabel('magnetic field [mT]');
