function [ok,data] = test(opt,olddata)

Exp.Sequence = 'MimsENDOR';
Exp.Field = 325;
Exp.tau = 0.1;
Exp.Range = larmorfrq('1H',Exp.Field) + [-1 1]*10;

Sys.Nucs = '1H';
Sys.A = [1 1 5];
Sys.lwEndor = 0.1;

Opt.Verbosity = 0;
[x,y] = saffron(Sys,Exp,Opt);
y = y/max(abs(y));

if (opt.Display)
  if ~isempty(olddata)
    subplot(3,1,[1 2]);
    plot(x,y,'r',x,olddata.y,'b');
    axis tight
    legend('new','old');
    title(mfilename);
    subplot(3,1,3);
    plot(x,olddata.y-y);
    axis tight
    xlabel('time (µs)');
  end
end

data.y = y;
data.x = x;

if isempty(olddata)
  ok = [];
else
  ok = areequal(olddata.y,y,1e-10,'rel');
end
