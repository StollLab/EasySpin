function ok = test()

% Intensities with higher-order Zeeman terms are the same with sparse and full
% matrices when the Hamiltonian is rediagonalized at the resonance fields
% (forced with Opt.RediagLimit = 1), including thermal populations

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A = 30;
Sys.Ham110 = -sqrt(3)*2*13.9962;  % (lB,lS,l) = (1,1,0), MHz/mT; equivalent to g = 2
Exp.mwFreq = 9.5;
Exp.Range = [300 380];
Exp.Temperature = 1;
Opt.RediagLimit = 1;

Opt.Sparse = false;
[Bfull,Ifull] = resfields(Sys,Exp,Opt);
Opt.Sparse = true;
[Bsparse,Isparse] = resfields(Sys,Exp,Opt);

ok = numel(Bfull)==2 && areequal(Bfull,Bsparse,1e-6,'abs') && ...
  areequal(Ifull,Isparse,1e-6,'rel');
