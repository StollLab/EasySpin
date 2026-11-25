% Orientation-selective two-pulse ESEEM (saffron)
%==========================================================================
% This example shows how saffron treats orientation selection: give an
% anisotropic g tensor (Sys.g), the microwave frequency (Exp.mwFreq) and
% the pulse lengths (Exp.tp). Opt.IdealPulses treats the pulses as ideal,
% with an excitation window of width 1/tp (here 200 MHz), which is faster
% than simulating the finite pulses. The two-pulse ESEEM spectrum is
% computed at a series of magnetic fields across the EPR spectrum.

clear, clf, clc

Exp.Sequence = '2pESEEM';
Exp.mwFreq = 9.5; % GHz
Exp.dt = 0.015; % µs
Exp.nPoints = 1024;
Exp.tau = 0.001; % µs
Exp.tp = 0.005; % µs, pulse lengths, excitation window 1/tp = 200 MHz

Sys.Nucs = '1H';
Sys.A = [1 6]; % MHz
Sys.g = [2 2.2];

Opt.GridSize = 91;
Opt.IdealPulses = true; % ideal pulses with an excitation window

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
