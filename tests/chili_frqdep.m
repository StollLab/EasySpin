function [err,data] = test(opt,olddata)

%=======================================================
% Frequency dependence of slow-motion spectra
%=======================================================
Sys = struct('g',[2.008 2.0061 2.0027],'Nucs','14N','A',[16 16 86]);
Sys.lw = 0.1;
Sys.tcorr = 3e-9;
Exp = struct('mwFreq',9.8);

mw = [3 9 15 35 95 140]; % GHz

for q = 1:numel(mw)
  Exp.mwFreq = mw(q);
  [x(q,:),y_] = chili(Sys,Exp);
  y(q,:) = y_/max(y_);
end

data.x = x;
data.y = y;

if isempty(olddata)
  err = [];
  return
end

ok = areequal(y,olddata.y,1e-4,'rel');
err = ~ok;

if opt.Display
  for k = 1:6
    subplot(3,2,k);
    plot(x(k,:),data.y(k,:),x(k,:),olddata.y(k,:),'.');
    axis tight
  end
end
