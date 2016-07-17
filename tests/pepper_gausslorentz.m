function [err,data] = test(opt,olddata)

% Test whether 
lwG = 1; lwL = 1;
Sys.S = 1/2;
Sys.g = 2;
Exp.mwFreq = 9.5;

for Harmonic = 0:2
  Exp.Harmonic = Harmonic;
  Sys.lw = 1;     [x,y] = pepper(Sys,Exp);
  Sys.lw = [1 0]; [x,y] = pepper(Sys,Exp);
  Sys.lw = [1 1]; [x,y] = pepper(Sys,Exp);
  Sys.lw = [0 1]; [x,y] = pepper(Sys,Exp);
end

err = zeros(1,12);
data = [];
