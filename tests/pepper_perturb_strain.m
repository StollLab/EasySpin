function ok = test(opt)

% Strain-broadened spectra with Opt.Method='perturb': agreement with 'matrix'
% where perturbation theory is accurate, and per-isotopologue strains

Optm.Method = 'matrix';
Optp.Method = 'perturb';

% Correlated g strain, field and frequency sweep
Sys = struct('g',[2.0 2.1 2.2],'StrainPars',{{'g(1)','g(3)'}},'StrainFWHM',[0.01 0.02],'StrainCorr',0.5);
Sys.lwpp = 0.3; % mT
Exp.mwFreq = 9.5;
Exp.Range = [300 350];
[B,y0] = pepper(Sys,Exp,Optm);
[~,y1] = pepper(Sys,Exp,Optp);
ok(1) = areequal(y0,y1,1e-4,'rel');
Sys.lwpp = 5; % MHz
ExpF.Field = 320;
ExpF.mwRange = [8.8 9.3];
[~,z0] = pepper(Sys,ExpF,Optm);
[~,z1] = pepper(Sys,ExpF,Optp);
ok(2) = areequal(z0,z1,1e-4,'rel');

% Correlated D/E strain for S = 5/2 at high field
Sys = struct('S',5/2,'D',[300 60],'StrainPars',{{'D(1)','D(2)'}},'StrainFWHM',[150 50],'StrainCorr',0.3);
Sys.lwpp = 0.5; % mT
Exp.mwFreq = 95;
Exp.Range = [3330 3460];
[B,y0] = pepper(Sys,Exp,Optm);
[~,y1] = pepper(Sys,Exp,Optp);
ok(3) = areequal(y0,y1,1e-3,'rel');

% Natural Cu with g/A strain: sum of the 63Cu and 65Cu isotopologues
Sys = struct('g',[2.05 2.25],'Nucs','Cu','A',[50 450],'lwpp',0.5, ...
  'StrainPars',{{'g(2)','A(2)'}},'StrainFWHM',[0.02 40],'StrainCorr',-0.5);
Exp.mwFreq = 9.5;
Exp.Range = [260 360];
[B,ynat] = pepper(Sys,Exp,Optp);
isos = {'63Cu','65Cu'};
abund = nucabund(isos);
[~,gn] = nucdata(isos);
yiso = 0;
for k = 1:2
  Sysk = Sys;
  Sysk.Nucs = isos{k};
  Sysk.A = Sys.A*gn(k)/gn(1);
  Sysk.StrainFWHM(2) = Sys.StrainFWHM(2)*gn(k)/gn(1);
  yiso = yiso + abund(k)*pepper(Sysk,Exp,Optp);
end
ok(4) = areequal(ynat,yiso,1e-6,'rel');

if opt.Display
  plot(B,ynat,B,yiso);
  legend('natural Cu','63Cu + 65Cu');
end
