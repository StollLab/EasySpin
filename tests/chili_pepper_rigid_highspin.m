function ok = test()

% Rigid-limit chili matches pepper in integrated intensity for S>1/2 and
% for two electrons, with and without temperature (frequency sweep)

Exp.Field = 339;  % mT
Exp.mwRange = [9.3 9.7];  % GHz
Exp.nPoints = 4096;
Exp.Harmonic = 0;

Opt.LLMK = [40 0 0 0];

Sys1.S = 1;
Sys1.g = 2;
Sys1.D = 100;  % MHz

Sys2.S = [1/2 1/2];
Sys2.g = [2 2.004; 2 2.004];  % axial, since chili needs an anisotropic system
Sys2.ee = 20;  % MHz

Systems = {Sys1,Sys2};
Temperatures = [NaN 300];  % K
for iSys = 1:numel(Systems)
  Sys = Systems{iSys};
  Sys.tcorr = 1e-5;  % s, rigid limit
  Sys.lw = [5 0.2];  % MHz
  for iT = 1:numel(Temperatures)
    Exp.Temperature = Temperatures(iT);
    [~,yp] = pepper(Sys,Exp);
    [~,yc] = chili(Sys,Exp,Opt);
    ok(iSys,iT) = areequal(sum(yc)/sum(yp),1,0.02,'abs');
  end
end
