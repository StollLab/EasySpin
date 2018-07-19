function [err,data] = test(opt,olddata)

%=======================================================
% lambda200 potential simulation
%=======================================================

Nx.Nucs = '14N';
Nx.g = [2.009 2.006 2.002];
Nx.A = mt2mhz([5 5.5 33]/10);
Nx.logDiff = 7;
Nx.lwpp = [0 0.0];

Exp.mwFreq = 9.54445;
Exp.CenterSweep = [340 12];
Exp.nPoints = 512;

Opt.LLKM = [10 5 2 2];
Opt.Verbosity = 0;
Opt.nKnots = 1;
Nx.Potential = [2 0 0 +1];

[x,y] = chili(Nx,Exp,Opt);
y = y.'/max(abs(y));

data.x = x;
data.y = y;

if ~isempty(olddata)
  if opt.Display
    subplot(4,1,[1,2,3]);
    plot(x,data.y,'r',x,olddata.y,'b');
    legend('new','old');
    subplot(4,1,4);
    plot(x,data.y-olddata.y);
  end
  ok = areequal(y,olddata.y,1e-4);
  err = ~ok;
else
  err = [];
end
