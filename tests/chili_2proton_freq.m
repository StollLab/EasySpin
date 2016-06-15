function [err,data] = test(opt,olddata)

%==========================================================================
% Frequency-swept slow-motion spectrum of a S=1/2 with hyperfine
%   splitting from 2 protons
%==========================================================================

Sys.S = 1/2;
Sys.g = [2.01 2.007 2.005];
Sys.A = [25 25];
Sys.Nucs = '1H,1H';
Sys.tcorr = 10e-9;

Exp.Field = 350;
Exp.mwRange = [9.65 10.2];

Opt.LiouvMethod = 'general';
[x,y] = chili(Sys,Exp,Opt);

data.x = x;
data.y = y;

% Check for consistency
if ~isempty(olddata)
  ok = areequal(y,olddata.y,1e-1);
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
