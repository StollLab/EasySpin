function [ok,data] = test(opt,olddata)

% Approximate field-sweep doublet, fast method

Sys.S = 1/2;
Sys.g = [2.01 2.005 2.002];
Sys.Nucs = '14N';
Sys.A = [20 20 100];
Sys.tcorr = 10e-9;

Exp.mwFreq = 9.5;
Exp.Range = [332 346];
Exp.Harmonic = 0;

Opt.LiouvMethod = 'fast';

[x,y] = chili(Sys,Exp,Opt);

data.x = x;
data.y = y;
 
% Check for consistency
if ~isempty(olddata)
  ok = areequal(y/max(y),olddata.y/max(olddata.y),1e-5,'abs');
else
  ok = [];
end

if opt.Display
  if ~isempty(olddata)
    subplot(4,1,1:3);
    plot(olddata.x,olddata.y/max(olddata.y),data.x,data.y/max(data.y));
    legend('old','new');
    legend boxoff
    subplot(4,1,4);
    plot(data.x,data.y/max(data.y)-olddata.y/max(olddata.y),'r');
    legend('new - old');
    legend boxoff
  end
end
