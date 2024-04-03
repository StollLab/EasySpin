function ok = test()

% Test whether Gaussian, Lorentzian, and Gaussian+Lorentzian convolutional
% broadenings work.

Sys.S = 1/2;
Sys.g = 2;
Exp.mwFreq = 9.5;
Exp.Range = [335 344];

lwG = 1;
lwL = 1;
for Harmonic = 0:2
  Exp.Harmonic = Harmonic;
  Sys.lw = lwG;       [B,spc] = pepper(Sys,Exp);
  Sys.lw = [lwG 0];   [B,spc] = pepper(Sys,Exp);
  Sys.lw = [lwG lwL]; [B,spc] = pepper(Sys,Exp);
  Sys.lw = [0   lwL]; [B,spc] = pepper(Sys,Exp);
end

ok = true;
