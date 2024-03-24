% Orientation-selective three-pulse ESEEM (saffron)
%==========================================================================
% This examples shows one way EasySpin can treat orientation selection.
% Just give a non-isotropic g tensor (Sys.g), a mw frequency (Exp.mwFreq)
% and an excitation width (Exp.ExcitationWidth).

clear, clf, clc

Exp.Sequence = '2pESEEM';
Exp.mwFreq = 9.5; % GHz
Exp.dt = 0.015; % mus
Exp.nPoints = 1024;
Exp.ExciteWidth = 200; % MHz
Exp.tau = 0.001; % mus

Sys.Nucs = '1H';
Sys.A = [1 6]; % MHz
Sys.g = [2 2.2];

Opt.GridSize = 91;

Fields = 300:2:350;

for k = 1:numel(Fields)
  Exp.Field = Fields(k);
  [t,td,info] = saffron(Sys,Exp,Opt);
  spc(k,:) = real(info.fd);
end

spc = spc/max(abs(spc(:)));
stackplot(info.f,spc,{'maxabs',4},Fields,compose('%1.0f mT',Fields));
xlim([0 max(info.f)]);
xlabel('\nu (MHz)')