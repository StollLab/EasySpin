function ok = test()

% Spectrum at negative fields is the mirror image of the spectrum at
% positive fields, for perturbation theory and the eigenfield method

Sys.g = [2 2.1 2.2];
Sys.lw = 1;
Exp.mwFreq = 9.5;
Exp.Harmonic = 0;

ok = true;
for method = {'perturb','eig'}
  Opt.Method = method{1};
  Exp.Range = [290 350];
  [~,yp] = pepper(Sys,Exp,Opt);
  Exp.Range = [-350 -290];
  [~,yn] = pepper(Sys,Exp,Opt);
  yr = fliplr(yp);
  % the mirrored spectrum is offset by at most one point
  ok = ok && (areequal(yn,yr,1e-8,'rel') || areequal(yn(1:end-1),yr(2:end),1e-8,'rel'));
end
