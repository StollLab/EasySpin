function ok = test(opt)

% Susceptibility of a simple coupled spin dimer

Sys.S = [1/2 1/2];
Sys.g = [2 2];
Sys.ee = -2*-4*30e3; % MHz

Exp.Temperature = 1:100; % K
Exp.Field = 1000; % mT

Opt.Output = 'chimol';
Opt.Method = 'operator';
chi1 = curry(Sys,Exp,Opt);
Opt.Method = 'partitionfunction';
chi2 = curry(Sys,Exp,Opt);

ok = areequal(chi1,chi2,1e-3,'rel');

if opt.Display
  plot(Exp.Temperature,chi1,Exp.Temperature,chi2);
  legend('operator','partitionfunction');
end
