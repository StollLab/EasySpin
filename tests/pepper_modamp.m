function [err,data] = test(opt,olddata)

Sys.lwpp = 0.5;
Sys.g = [2 2.02];
Exp.mwFreq = 9.7;
Exp.CenterSweep = [346 20];
Exp.ModAmp = 5;
Opt = struct('Verbosity',opt.Verbosity);
[x,y] = pepper(Sys,Exp,Opt);

data.x = x;
data.y = y;

if ~isempty(olddata)
  delta = 1e-6*max(abs(olddata.y));
  err = ~areequal(y,olddata.y,delta);
else
  err = [];
end

if (opt.Display)
  if ~isempty(olddata)
    subplot(4,1,[1 2 3]);
    plot(olddata.x,olddata.y,'b',data.x,data.y,'r');
    legend('old','new');
    subplot(4,1,4);
    plot(data.x,olddata.y-data.y);
  else
    plot(x,y);
  end
  xlabel('magnetic field [mT]');
  axis tight;
end
