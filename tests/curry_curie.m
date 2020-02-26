function ok = test(opt)

% Compare molar susceptibility calculated with curry() against Curie law

S = 1/2;
g = 2;
T = 1:100; % K
B = 0; % T
chimol_curie = S*(S+1)*g^2*mu0*bmagn^2*avogadro/3/boltzm./T;

Sys.S = S;
Sys.g = g;
Exp.Temperature = T;
Exp.Field = B;
Opt.Output = 'chimol';
chimol_curry = curry(Sys,Exp,Opt);

if opt.Display
  plot(T,chimol_curie,T,chimol_curry);
  xlabel('temperature (K)');
  ylabel('molar magnetic susceptibility  (K m^3 mol^{-1})');
  legend('Curie law','curry');
  legend boxoff
end

ok = areequal(chimol_curie,chimol_curry,1e-10,'abs');

