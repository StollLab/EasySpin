function [err,data] = test(opt,olddata)

%=======================================================
% Frequency-sweep doublet, Freed solver
%=======================================================

Sys.g = [2.01 2];
Sys.Nucs = '14N';
Sys.A = [20 20 100];
Sys.tcorr = 10e-9;

Exp.Field = 350;
Exp.mwRange = [9.6 10];
Opt.LiouvMethod = 'Freed';

[x,y] = chili(Sys,Exp,Opt);

data.x = x;
data.y = y;

% Check for consistency
if ~isempty(olddata)
  ok = areequal(y,olddata.y,1e-2*max(y));
  err = ~ok;
else
  err = [];
end

if opt.Display
  if ~isempty(olddata)
    subplot(4,1,1:3);
    plot(data.x,data.y,'r',olddata.x,olddata.y,'g');
    legend('new','old');
    legend boxoff
    subplot(4,1,4);
    plot(data.x,data.y-olddata.y,'r');
    legend('new - old');
    legend boxoff
  end
end
