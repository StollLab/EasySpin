function [err,data] = test(opt,olddata)

% Susceptibility of a simple coupled spin dimer

Sys.S = [1/2 1/2];
Sys.g = [2 2];
Sys.ee = -2*-4*30e3; % MHz

Exp.Temperature = 1:100; % K
Exp.Field = 1000; % mT

Opt.Output = 'chimol';
chi_SI = curry(Sys,Exp,Opt);

data.chi_SI = chi_SI;

if opt.Display
  T = Exp.Temperature;
  plot(T,data.chi_SI,T,olddata.chi_SI);
  legend('new','old');
  xlabel('temperature (K)');
  ylabel('magnetic susceptibility  (K cm^3 mol^{-1})');
end

if ~isempty(olddata)
  ok = areequal(data.chi_SI,olddata.chi_SI,max(data.chi_SI)*1e-3,'abs');
  err = ~ok;
else
  err = [];
end
