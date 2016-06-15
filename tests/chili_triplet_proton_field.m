function [err,data] = test(opt,olddata)

%==========================================================================
% Field-swept slow-motion spectrum of a triplet with hyperfine
%   splitting from one proton
%==========================================================================

Sys.S = 1;
Sys.D = 100;
Sys.A = 10;
Sys.Nucs = '1H';
Sys.tcorr = 10e-9;

Exp.mwFreq = 9.5;
%Exp.CenterField = 1e3*planck*Exp.mwFreq*1e9/(bmagn*mean(Sys.g));
Exp.Range = [320 360];
Exp.Harmonic = 0;

[x,y] = chili(Sys,Exp);

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
