function [ok,data] = test(opt,olddata)

% Regression test: Two coupled electrons, two D strains
%------------------------------------------------------

Sys.S = [1 1];
Sys.g = [2 2.5];
Sys.D = [400 700];
Sys.ee = 1;
Sys.lwpp = 1;

Exp.mwFreq = 9.5;
Exp.Range = [200 400];

Ds1 = [30 10];
Ds2 = [100 20];
Sys.DStrain = [Ds1; Ds2];
[x,spc] = pepper(Sys,Exp);

if (opt.Display)
  if isempty(olddata)
    plot(x,spc);
  else
    h = plot(x,spc,'k',x,olddata.spc,'r');
    legend('old','new');
    %set(h(1),'LineWidth',2);
  end
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: S1=S2=1 with two D strains');
end

data.spc = spc;

if ~isempty(olddata)
  ok = areequal(spc,olddata.spc,1e-3,'rel');
else
  ok = [];
end
