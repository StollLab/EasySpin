function ok = test(opt)

% Check that a Lorentzian narrower than two increments does not prevent
% computing a harmonic when a Gaussian width is present

Sys.g = [2 2.05 2.1];
Exp.mwFreq = 9.5;
Exp.Range = [315 345];
Exp.nPoints = 2000;
Exp.Harmonic = 1;

dx = diff(Exp.Range)/(Exp.nPoints-1);

Sys.lw = [1 0];
[x,spc0] = pepper(Sys,Exp);
Sys.lw = [1 dx/4];
[x,spc1] = pepper(Sys,Exp);

ok = areequal(spc1,spc0,1e-2,'rel');

if opt.Display
  plot(x,spc0,x,spc1,'--');
  legend('lw = [1 0]','lw = [1 dx/4]');
end

end
