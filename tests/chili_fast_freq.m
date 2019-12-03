function [err,data] = test(opt,olddata)

%=======================================================
% Frequency-sweep doublet, fast method
%=======================================================

Sys.g = [2.01 2];
Sys.Nucs = '14N';
Sys.A = [20 20 100];
Sys.tcorr = 10e-9;

Exp.Field = 350;
Exp.mwRange = [9.6 10];
Opt.LiouvMethod = 'general';
Opt.Solver = 'L';

[x,y] = chili(Sys,Exp,Opt);

data.x = x;
data.y = y;

% Check for consistency
if ~isempty(olddata)
  ok = areequal(y/max(y),olddata.y/max(olddata.y),1e-2,'abs');
  err = ~ok;
else
  err = [];
end

if opt.Display
  if ~isempty(olddata)
    subplot(4,1,1:3);
    y_rescaled = data.y/max(data.y)*max(olddata.y);
    plot(olddata.x,olddata.y,data.x,data.y,data.x,y_rescaled);
    legend('old','new','rescaled new');
    legend boxoff
    subplot(4,1,4);
    plot(data.x,data.y-olddata.y,'r');
    legend('new - old');
    legend boxoff
  end
end
