function ok = test()

% Negative magnetic field is rejected

Sys.g = [2 2.1 2.2];
Sys.lw = 1;
Sys.tcorr = 1e-9;
Exp.Field = -340;
Exp.mwCenterSweep = [9.5 1];
Par.Model = 'diffusion';
Par.dtSpin = 1e-9;
Par.nSteps = 10;
Par.nTraj = 1;

try
  cardamom(Sys,Exp,Par);
  ok = false;
catch err
  ok = contains(err.message,'Exp.Field cannot be negative');
end
