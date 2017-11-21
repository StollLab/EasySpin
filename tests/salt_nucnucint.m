function [err,data] = test(opt,olddata)

%===============================================================================
% Make sure nuclear-nuclear splitting appears in ENDOR spectrum
%===============================================================================

Sys.Nucs = '1H,1H';
Sys.A = [5 10]; % MHz
Sys.nn = 0.2;
Sys.lwEndor = 0.1; % MHz

Exp.Field = 1200; % mT
Exp.Range = [40 60]; % MHz
Exp.nPoints = 2000;

Opt.Method = 'matrix';

[nu,spc] = salt(Sys,Exp,Opt);

if opt.Display
  plot(olddata.nu,olddata.spc,nu,spc);
  legend('old','new');
  legend boxoff
end

data.nu = nu;
data.spc = spc;

if ~isempty(olddata)
  err = ~areequal(olddata.spc,spc,1e-5);
else
  err = [];
end

