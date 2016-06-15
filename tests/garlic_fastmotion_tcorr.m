function [err,data] = test(opt,olddata)

%=======================================================
% Fast motional regime: specifying the correlation time
%=======================================================

Spins.g = [2.0088 2.0061 2.0027];
Spins.Nucs = '14N';
Spins.A = mt2mhz([5.8 5.8 30.8]/10); % G -> MHz

Spins.tcorr = 2e-9;

Exp.mwFreq = 9.7;
Exp.Range = [341 350];
[x,y] = garlic(Spins,Exp);

data.x = x;
data.y = y;

if ~isempty(olddata)
  ok = areequal(x,olddata.x) & areequal(y,olddata.y,1e-6*max(abs(y)));
  err = ~ok;
else
  err = [];
end

if (opt.Display)
  subplot(2,1,1);
  plot(olddata.x,olddata.y,'k',x,data.y,'r');
  xlabel('magnetic field [mT]'); axis tight
  subplot(2,1,2);
  plot(x,data.y-olddata.y,'r');
  pause;
end
