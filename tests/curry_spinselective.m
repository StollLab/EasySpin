function [ok,data] = test(opt,olddata)

% Test spin-selective magnetic moments

Sys.S = [7/2 1/2];
Sys.ee = -2*5*30e3; % MHz

T = 4; % K
B = 0:10:3000; % mT

Exp.Temperature = T;
Exp.Field = B;
Opt.Output = 'mu';

Opt.Spins = 1;
mu1 = curry(Sys,Exp,Opt);
Opt.Spins = 2;
mu2 = curry(Sys,Exp,Opt);

data.mu1 = mu1;
data.mu2 = mu2;

if ~isempty(olddata)
  thr = 1e-4*max(mu1);
  ok(1) = areequal(data.mu1,olddata.mu1,thr,'abs');
  ok(2) = areequal(data.mu2,olddata.mu2,thr,'abs');
else
  ok = [];
end
