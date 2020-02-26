function ok = chili_magnetictilt(opt)

% Invariance of spectrum when tilting A w.r.t. isotropic g

Sys.g = 2.0088;
Sys.Nucs = '14N';
Sys.A = [17 17 90];
Sys.Diff = 6e7;
Exp.mwFreq = 9.5;
Exp.Range = [330.0 345.0];
Opt = struct;

Sys.AFrame = [0,0,0];
[x,y0] = chili(Sys,Exp,Opt);
Sys.AFrame = [0,30,0]*pi/180;
[x,y1] = chili(Sys,Exp,Opt);

y0 = y0./max(y0);
y1 = y1./max(y1);

if opt.Display
  plot(x,y0,'-r',x,y1,'-b');
  legend('no tilt','tilt');
end

ok = areequal(y0,y1,2e-3,'abs');
