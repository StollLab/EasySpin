function [err,data] = test(opt,olddata)

%=======================================================
% Explicit field-sweep triplet, general solver
%=======================================================

Sys.S = 1;
Sys.g = [2.01 2.005 2.002];
Sys.D = 100;
Sys.tcorr = 10e-9;

Exp.mwFreq = 9.5;
Exp.Range = [332 345];
Exp.Harmonic = 0;

Opt.LiouvMethod = 'general';
Opt.FieldSweepMethod = 'explicit';
Opt.LLMK = [6 0 0 2];
Opt.pSmin = 1;

[x,y] = chili(Sys,Exp,Opt);

data.x = x;
data.y = y;

% Check for consistency
if ~isempty(olddata)
  ok = areequal(y/max(y),olddata.y/max(olddata.y),1e-1,'abs');
  err = ~ok;
else
  err = [];
end

if opt.Display
  if ~isempty(olddata)
    subplot(4,1,1:3);
    plot(data.x,data.y/max(data.y),'r',olddata.x,olddata.y/max(olddata.y),'g');
    legend('new','old');
    legend boxoff
    subplot(4,1,4);
    plot(data.x,data.y/max(data.y)-olddata.y/max(olddata.y),'r');
    legend('new - old');
    legend boxoff
  end
end
