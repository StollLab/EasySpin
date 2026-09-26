function ok = test()

% Sys.initState is not supported and must give an error

Sys.S = 1;
Sys.D = 300;  % MHz
Sys.initState = {[1 0 0],'zerofield'};

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [250 420];  % mT
Exp.SampleFrame = [0 0 0];

try
  resfields_perturb(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'Sys.initState is not supported');
end
