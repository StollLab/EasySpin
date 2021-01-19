function [ok,data] = test(opt,olddata)

Exp.Sequence = '2pESEEM';
Exp.mwFreq = 9.5;
Exp.dt = 0.015;
Exp.nPoints = 128;
Exp.ExciteWidth = 200;
Exp.tau = 0.001;
Exp.Field = 350;

Sys.Nucs = '1H';
Sys.A = [1 6];
Sys.g = [2 2.2];

Opt.GridSize = 91;

[x,y] = saffron(Sys,Exp,Opt);

data.x = x;
data.y = y;

if (opt.Display)
  if ~isempty(olddata)
    subplot(3,1,[1 2]);
    plot(x',real(y)','r',x',real(olddata.y)','b');
    axis tight
    legend('new','old');
    title(mfilename);
    subplot(3,1,3);
    plot(x',real(olddata.y-y)');
    axis tight
    xlabel('time [us]');
  end
end

if ~isempty(olddata)
  ok = areequal(y,olddata.y,1e-4,'abs');
else
  ok = [];
end

