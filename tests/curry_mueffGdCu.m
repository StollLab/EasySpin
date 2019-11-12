function [err,data] = test(opt,olddata)

% Test calculation of effective magnetic moment for a Gd(III)-Cu(II) complex.
% This reproduces Fig. 1 from Hatscher et al, Pure Appl. Chem. 77, 497 (2005).

Sys.S = [7/2 1/2];
Sys.ee = -2*5*30e3; % MHz

T = 0:0.1:125; % K

Exp.Temperature = T;

Opt.Output = 'mueff';
Opt.deltaB = 0.00001;

Exp.Field = 10; % mT
mueff1 = curry(Sys,Exp,Opt);
Exp.Field = 100; % mT
mueff2 = curry(Sys,Exp,Opt);
Exp.Field = 1000; % mT
mueff3 = curry(Sys,Exp,Opt);

data.mueff1 = mueff1;
data.mueff2 = mueff2;
data.mueff3 = mueff3;

if opt.Display
  plot(T,mueff1,T,mueff2,T,mueff3,...
    T,olddata.mueff1,'.',T,olddata.mueff2,'.',T,olddata.mueff3,'.');
  axis tight
  set(gca,'XTick',0:25:125);
  ylim([0 9]);
  
  xlabel('temperature (K)');
  ylabel('\mu_{eff}');
  legend('B_0 = 0.01 T','B_0 = 0.1 T','B_0 = 1 T','location','best');
  legend boxoff
end

if ~isempty(olddata)
  thr = 1e-4*max(data.mueff1);
  ok(1) = areequal(olddata.mueff1,data.mueff1,thr,'abs');
  ok(2) = areequal(olddata.mueff2,data.mueff2,thr,'abs');
  ok(3) = areequal(olddata.mueff3,data.mueff3,thr,'abs');
  err = any(~ok);
else
  err = [];
end
