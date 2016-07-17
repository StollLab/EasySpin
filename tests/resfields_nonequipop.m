function [err,data] = test(opt,olddata)

% Test correct sign of emission/absorption for non-Boltzmann populations
% Low-field line must be in emission, the high-field line in absorption
Sys.S = 1;
Sys.g = [2 2 2];
D = 0.06; E = 0;
Sys.D = 100*clight/1e6*[-D/3 + E, -D/3 - E, 2*D/3];

Exp = struct('mwFreq',9.67739,'nPoints',426,'Range',[260 430]);
Exp.Harmonic = 0;
Exp.Temperature = [1 0 0];
Exp.CrystalOrientation = [0 0 0];

[B,A] = resfields(Sys,Exp);

[B,idx] = sort(B);
A = A(idx);

data = [];

if numel(B)==2
  err = (A(1)<=0) | (A(2)>=0);
else
  err = 1;
end
