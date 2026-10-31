function [ok,data] = test(opt,olddata)

% Regression test: 3-pulse ESEEM with T1 and T2 decay (fast mode)

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2];  % MHz
Sys.T1 = 2;  % µs
Sys.T2 = 1;  % µs

Opt.GridSize = 20;

Exp.Sequence = '3pESEEM';
Exp.dt = 0.01;  % µs
Exp.tau = 0.1;  % µs
Exp.nPoints = 191;
Exp.Field = 350;  % mT

[T, V] = saffron(Sys,Exp,Opt);
V = real(V);

if opt.Display && ~isempty(olddata)
  plot(T,V,olddata.T,olddata.V);
  legend('current','reference');
  xlabel('\tau + T (µs)');
  ylabel('V');
  grid on
end

data.T = T;
data.V = V;

if ~isempty(olddata)
  ok = areequal(T,olddata.T,1e-5,'rel') && areequal(V,olddata.V,1e-5,'rel');
else
  ok = [];
end
