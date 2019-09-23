function [err,data] = test(opt,olddata)

% Compare powder simulation of isotropic system to crystal simulation

Sys.S = [1/2 1/2];
Sys.g = [2 2.15];
Sys.ee = 400; % MHz
Sys.lw = 30;

Exp.Field = 350;
Exp.Harmonic = 0;
Exp.mwRange = [9.2 11];
Exp.nPoints = 50000;

% (1) Powder averaging (zone projection)
Opt.Symmetry = 'Dinfh';
[x,y2] = pepper(Sys,Exp,Opt);

% (2) Powder averaging (triangle projection)
Opt.Symmetry = 'Ci'; 
[x,y1] = pepper(Sys,Exp,Opt);

% (3) Single-crystal (accumulation using template copy)
Exp.CrystalOrientation = [0 0 0];
[x,y0] = pepper(Sys,Exp,Opt);

if opt.Display
  plot(x,y0,x,y1,'r',x,y2,'g');
  legend('crystal','powder Ci','powder Dinfh');
  legend boxoff
end

thresh = 1e-2;
err = ~areequal(y0,y1,thresh,'rel') || ~areequal(y0,y2,thresh,'rel');
data = [];
