function [ok,data] = test(opt,olddata)

% Check whether resfreqs_matrix handles gStrain correctly.

Sys.S = 1/2;
Sys.g = 2;
Sys.gStrain = 0.01;

mw = 9.8; % GHz
Exp.Field = mw*1e9*planck/bmagn/Sys.g/1e-3;  % mT
Exp.mwRange = mw + [-1 1]*0.1;  % GHz

[B,spc] = pepper(Sys,Exp);
spc = spc/max(spc);

if opt.Display
  plot(B,spc,B,olddata.y);
  yline(0.5);
  xline(mw);
  g_up = Sys.g + Sys.gStrain/2;
  g_lo = Sys.g - Sys.gStrain/2;
  dmw_up = g_up*bmagn*Exp.Field*1e-3/planck/1e9;
  dmw_lo = g_lo*bmagn*Exp.Field*1e-3/planck/1e9;
  xline(dmw_up);
  xline(dmw_lo);
  legend('new','old');
end
data.y = spc;

if ~isempty(olddata)
  ok = areequal(spc,olddata.y,1e-3,'abs');
else
  ok = [];
end
