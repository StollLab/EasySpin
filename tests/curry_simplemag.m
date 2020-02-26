function [ok,data] = test(opt,olddata)

% Magnetization of a simple coupled spin dimer

Sys.S = [1/2 1/2];
Sys.g = [2 2];
Sys.ee = -2*-4*30e3; % MHz

Exp.Temperature = 1; % K
Exp.Field = linspace(0,17,20)*1e3; % mT

Opt.Output = 'muBM'; % identical for SI and CGS-emu
m = curry(Sys,Exp,Opt);

data.m = m;

if opt.Display
  B = Exp.Field/1e3;
  plot(B,olddata.m,B,data.m);
  legend('old','new');
  legend boxoff
  xlabel('magnetic field (T)');
  ylabel('magnetic moment');
end

if ~isempty(olddata)
  ok = areequal(data.m,olddata.m,1e-3,'abs');
else
  ok = [];
end
