function [err,data] = test(opt,olddata)

Exp.Sequence = 'MimsENDOR';
Exp.Field = 325;
Exp.tau = 0.1;
Exp.Range = larmorfrq('1H',Exp.Field) + [-1 1]*10;

Sys1.Nucs = '1H';
Sys1.A = [1 3];
Sys1.lwEndor = 0.1;
Sys2.Nucs = '1H';
Sys2.A = [7 9];
Sys2.lwEndor = 0.1;
Sys2.weight = 0.3;

Opt.Verbosity = 0;
[x,y] = saffron({Sys1,Sys2},Exp,Opt);
y = y/max(abs(y));

if (opt.Display)
  if ~isempty(olddata)
    subplot(3,1,[1 2]);
    plot(x,y,'r',x,olddata.y,'b');
    xlabel('time [us]');
    legend('new','old');
    title(mfilename)
    subplot(3,1,3);
    plot(x,olddata.y-y);
    title('old - new');
  end
end

data.y = y;
data.x = x;

if isempty(olddata)
  err = [];
else
  ok = areequal(olddata.y,y,1e-10,'rel');
  err = ~ok;
end
