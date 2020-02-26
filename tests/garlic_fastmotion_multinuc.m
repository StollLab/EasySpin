function [ok,data] = test(opt,olddata)

%=======================================================
% Fast motional regime: specifying the correlation time
%=======================================================

Spins.g = [2.0088 2.0061 2.0027];
Spins.Nucs = '14N,1H';
Spins.A = mt2mhz([5.8 5.8 30.8; -3 -3 10]/10); % G -> MHz

Spins.tcorr = 1e-9;

Exp = struct('mwFreq',9.7);
Exp.Range = [343 348];
[x,y] = garlic(Spins,Exp);

data.y = y;

if ~isempty(olddata)
  ok = areequal(y,olddata.y,1e-7,'rel');
else
  ok = [];
end

if opt.Display
  subplot(4,1,[1 2 3]);
  plot(x,olddata.y,'k',x,data.y,'r');
  xlabel('magnetic field [mT]'); axis tight
  subplot(4,1,4);
  plot(x,olddata.y-data.y);
end
