function [err,data] = test(opt,olddata)

% Regression test: correlated D strain

Sys.S = 1;
Sys.D = [800 80];
Sys.DStrain = [100 33];
Sys.lwpp = 1;

Exp.mwFreq = 9.5;
Exp.CenterSweep = [340 100];

Opt.nKnots = 15;

Sys.DStrainCorr =  0; [x,y0] = pepper(Sys,Exp,Opt);
Sys.DStrainCorr = +1; [x,yp] = pepper(Sys,Exp,Opt);
Sys.DStrainCorr = -1; [x,ym] = pepper(Sys,Exp,Opt);

if (opt.Display)
  if isempty(olddata)
    plot(x,y0,'b',x,yp,'r',x,ym,'g');
    legend('corr = 0','corr = +1','corr = -1');
    legend boxoff
  else
    subplot(3,1,1); plot(x,olddata.y0,'b',x,y0,'r');
    subplot(3,1,2); plot(x,olddata.yp,'b',x,yp,'r');
    subplot(3,1,3); plot(x,olddata.ym,'b',x,ym,'r');
    legend('old','new');
    legend boxoff;
  end
end

data.y0 = y0;
data.yp = yp;
data.ym = ym;

if ~isempty(olddata)
  m = max([y0 yp ym]);
  err = ~areequal(y0,olddata.y0,m*0.01);
  err = err || ~areequal(yp,olddata.yp,m*0.01);
  err = err || ~areequal(ym,olddata.ym,m*0.01);
else
  err = [];
end
