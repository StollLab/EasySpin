function [ok,data] = test(opt,olddata)

Sys.Nucs = '1H';
Sys.A = [2 9];

Exp.Sequence = '4pESEEM';
Exp.tau = 0.100;
Exp.dt = 0.010;
%Exp.T = 0.080;
Exp.Field = 350;

Opt.Verbosity = 0;
[x,y] = saffron(Sys,Exp,Opt);
y = y/max(abs(y));

if (opt.Display)
  xlabel('time [us]');
  if ~isempty(olddata)
    subplot(3,1,[1 2]);
    plot(x,y,'r',x,olddata.y,'b');
    legend('new','old');
    title(mfilename);
    subplot(3,1,3);
    plot(x,y-olddata.y,'r');
    title('residuals');
  end
end

data.y = y;
data.x = x;

if isempty(olddata)
  ok = [];
else
  ok = areequal(olddata.y,y,1e-10,'rel');
end
