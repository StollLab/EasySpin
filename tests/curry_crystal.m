function [ok,data] = test(opt,olddata)

% Check single-crystal calculation of magnetic susceptibility
% (reproduces Kahn, Molecular Magnetism, p. 19, Figure 2.4 top)

D = 5*30e3; % MHz
T = 0:0.1:50; % K

Sys.S = 1;
Sys.g = 2;
Sys.D = D;

Exp.Temperature = T;
Exp.Field = 0;

Opt.Units = 'CGS';
Opt.Output = 'chimol';

Exp.CrystalOrientation = [0 0 0];
chiz = curry(Sys,Exp,Opt);
Exp.CrystalOrientation = [0 pi/2 0];
chix = curry(Sys,Exp,Opt);

data.chix = chix;
data.chiz = chiz;

if opt.Display
  h = plot(T,data.chix,T,data.chiz,T,olddata.chix,T,olddata.chiz);
  set(h(3:4),'LineStyle','--');
  xlabel('{\itT} (K)');
  ylabel('{\it\chi}_u (cm^3 mol^{-1})');
  legend('\chi_x','\chi_z','\chi_x old','\chi_z old');
  legend boxoff
end

if ~isempty(olddata)
  thr = 1e-5;
  ok = areequal(olddata.chix,chix,thr,'abs') && areequal(olddata.chiz,chiz,thr,'abs');
else
  ok = [];
end
