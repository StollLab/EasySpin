function [err,data] = test(opt,olddata)

clear Sys Exp Opt
Exp.Sequence = '2pESEEM';
Exp.Field = 324.9;
Exp.dt = 0.050;
Exp.nPoints = 1001;
Exp.tau = 0.001;

Sys.Nucs = '14N';
nuI = larmorfrq(Sys.Nucs,Exp.Field);
Delta = +0;
Sys.A = (2+Delta)*nuI + [-1 2]*0*nuI;
Sys.Q = [4*0.1*nuI, 0.6];

Opt.Verbosity = 0;
[x,y] = saffron(Sys,Exp,Opt);
y = y/max(abs(y));

if (opt.Display)
  if ~isempty(olddata)
    subplot(3,1,[1 2]);
    plot(x,y,'r',x,olddata.y,'b');
    xlabel('frequency [MHz]');
    legend('new','old');
    subplot(3,1,3);
    plot(x,y-olddata.y);
    title('new - old');
  end
end

data.y = y;
data.x = x;

if isempty(olddata)
  err = [];
else
  ok = areequal(olddata.y,y);
  err = ~ok;
end
