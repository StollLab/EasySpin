function ok = test(opt)

% Verify projective approximation for various
% values of relative spectral broadening

Sys = struct('S',1/2,'g',[2 2.4 3]);
Exp = struct('mwFreq',9.5,'Range',[200 400],'Harmonic',0);
Opt = struct('Verbosity',opt.Verbosity,'nKnots',[19 3],'Method','perturb');
Exp.nPoints = 1e4;

lw = [1 3 10 30 100 300];
for k=1:numel(lw)
  Sys.HStrain = [1 1 1]*lw(k);
  [x,y(k,:)] = pepper(Sys,Exp,Opt);
end

if (opt.Display)
  plot(x,y);
  title('Various relative broadenings, anisotropic');
end

ok = true;
