function [ok,data] = test(opt,olddata)

% Regression test: correlated D strain

Sys.S = 1;
Sys.D = [800 80];
Sys.DStrain = [100 33];
Sys.lwpp = 1;

Exp.mwFreq = 9.5;
Exp.CenterSweep = [340 100];

Opt.GridSize = 15;

Sys.DStrainCorr =  0; [B,spc0] = pepper(Sys,Exp,Opt);
Sys.DStrainCorr = +1; [B,spcp] = pepper(Sys,Exp,Opt);
Sys.DStrainCorr = -1; [B,spcm] = pepper(Sys,Exp,Opt);

if opt.Display
  if isempty(olddata)
    plot(B,spc0,'b',B,spcp,'r',B,spcm,'g');
    legend('corr = 0','corr = +1','corr = -1');
    legend boxoff
  else
    subplot(3,1,1); plot(B,olddata.y0,'b',B,spc0,'r');
    title('corr = 0');
    subplot(3,1,2); plot(B,olddata.yp,'b',B,spcp,'r');
    title('corr = +1');
    subplot(3,1,3); plot(B,olddata.ym,'b',B,spcm,'r');
    title('corr = -1');
    legend('old','new');
    legend boxoff;
  end
end

data.y0 = spc0;
data.yp = spcp;
data.ym = spcm;

if ~isempty(olddata)
  m = max([spc0 spcp spcm]);
  ok = areequal(spc0,olddata.y0,m*0.01,'abs');
  ok = ok && areequal(spcp,olddata.yp,m*0.01,'abs');
  ok = ok && areequal(spcm,olddata.ym,m*0.01,'abs');
else
  ok = [];
end
