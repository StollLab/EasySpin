function [err,data] = test(opt,olddata)

Sys.lwpp = 0.5;
Exp.mwFreq = 9.7;
Exp.CenterSweep = [346 10];
Exp.ModAmp = 5;
Opt = struct('Verbosity',opt.Verbosity);
[x,y] = garlic(Sys,Exp,Opt);

data.x = x;
data.y = y;

if ~isempty(olddata)
  err = ~areequal(x,olddata.x) | ~areequal(y,olddata.y);
else
  err = [];
end

if (opt.Display)
  if ~isempty(olddata)
    plot(olddata.x,olddata.y,'b',data.x,data.y,'r');
    legend('old','new');
  else
    plot(x,y);
  end
  xlabel('magnetic field [mT]');
  axis tight;
end
