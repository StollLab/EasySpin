function [err,data] = test(opt,olddata)

% Test whether Gaussian, Lorentzian, and Gaussian+Lorentzian broadenings always
% give requested harmonic.

Sys.g = 2;
Exp.mwFreq = 9.5;

lwppG = 0.1; lwppL = 0.1;

isfirstharmonic = @(s)areequal(max(s),-min(s),1e-3,'abs');

% Gaussian only
Sys.lwpp = lwppG;
[B,spc1] = pepper(Sys,Exp);
err(1) = ~isfirstharmonic(spc1);

% Gaussian with zero Lorentzian
Sys.lwpp = [lwppG 0];
[B,spc2] = pepper(Sys,Exp);
err(2) = ~isfirstharmonic(spc2);

% Lorentzian with zero Gaussian
Sys.lwpp = [0 lwppL];
[B,spc3] = pepper(Sys,Exp);
err(3) = ~isfirstharmonic(spc3);

% Combined Gaussian and Lorentzian
Sys.lwpp = [lwppG lwppL];
[B,spc4] = pepper(Sys,Exp);
err(4) = ~isfirstharmonic(spc4);

data = [];
