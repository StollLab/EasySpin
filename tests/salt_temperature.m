function [err,data] = test(opt,olddata)

Sys.g = 1.75;
Sys.Nucs = '1H';
Sys.lwEndor = 0.055; % in MHz

Exp.mwFreq = 94.26;
Exp.Field = 3850;
Exp.Range = [161.5 166.5];
Exp.ExciteWidth = 222.7;
Exp.Temperature = 10;

Sys.A = -0.8;
[x,y1] = salt(Sys,Exp);

Sys.A = -Sys.A;
[x,y2] = salt(Sys,Exp);


if (opt.Display)
  subplot(3,1,[1 2]);
  h = plot(x,y1,'b',x,y2,'r');
  set(h(1),'LineWidth',2);
  legend('negative A','positive A');
  xlabel('frequency [MHz]');
  if ~isempty(olddata)
    subplot(3,1,3);
    plot(x,y1-olddata.y1,'b',x,y2-olddata.y2,'r');
    legend('y1 residuals','y2 residuals');
  end
end

data.y1 = y1;
data.y2 = y2;

if isempty(olddata)
  err = [];
else
  ok = ...
    areequal(y1,olddata.y1,1e-6,'rel') & ...
    areequal(y2,olddata.y2,1e-6,'rel');
  err = ~ok;
end
