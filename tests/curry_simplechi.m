function [err,data] = test(opt,olddata)

% Susceptibility of a simple coupled spin dimer

Sys.S = [1/2 1/2];
Sys.g = [2 2];
Sys.ee = -2*-4*30e3;

Exp.Temperature = 1:100;
Exp.Field = 1000;

[ignore,chi_SI] = curry(Sys,Exp);

data.chi_SI = chi_SI;

if opt.Display
  T = Exp.Temperature;
  subplot(1,2,1);
  plot(T,data.chi_SI,T,olddata.chi_SI);
  legend('new','old');
  xlabel('temperature (K)');
  ylabel('magnetic susceptibility  (K cm^3 mol^{-1})');
  subplot(1,2,2);
  c = 4*pi*1e-6;
  plot(T,data.chi_SI/c,T,olddata.chi_SI/c);
  legend('new','old');
  xlabel('temperature (K)');
  ylabel('magnetic susceptibility  (K cm^3 mol^{-1})');
end

if ~isempty(olddata)
  ok = areequal(data.chi_SI,olddata.chi_SI,1e-12);
  err = ~ok;
else
  err = [];
end
