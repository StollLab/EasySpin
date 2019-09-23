function [err,data] = test(opt,olddata)

%=======================================================
% Frequency dependence of slow-motion spectra
%=======================================================
Sys = struct('g',[2.008 2.0061 2.0027],'Nucs','14N','A',[16 16 86]);
Sys.lw = 0.1;
Sys.tcorr = 1e-8;
Exp = struct('mwFreq',9.8);

mw = [3 9 15 35 95 140]; % GHz

for q = 1:numel(mw)
  Exp.mwFreq = mw(q);
  [x(q,:),y(q,:)] = chili(Sys,Exp);
  y(q,:) = y(q,:)/max(y(q,:));
end

data.x = x;
data.y = y;

if ~isempty(olddata)
  if opt.Display
    xx = 1:size(x,2);
    plot(xx,data.y,xx,olddata.y,'.');
    legend('new','old');
  end
  ok = areequal(y,olddata.y,1e-4,'abs');
  err = ~ok;
else
  err = [];
end
