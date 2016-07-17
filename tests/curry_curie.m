function [err,data] = test(opt,olddata)

% Compare curry susceptibility against Curie law

S = 1/2;
g = 2;
T = 1:100;
B = 0;
chi_curie = S*(S+1)*g^2*mu0*bmagn^2*avogadro/3/boltzm./T;

Sys.S = S;
Sys.g = g;
Exp.Temperature = T;
Exp.Field = B;
[muz,chi_curry] = curry(Sys,Exp);

if opt.Display
  plot(T,chi_curie,T,chi_curry);
  legend('Curie','curry');
  xlabel('temperature (K)');
  ylabel('magnetic susceptibility  (K m^3 mol^{-1})');
end

ok = areequal(chi_curie,chi_curry,1e-12);
err = ~ok;

data = [];
