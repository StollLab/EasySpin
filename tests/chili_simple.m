function [err,data] = test(opt,olddata)

%=======================================================
% Simple isotropic simulations
%=======================================================
Sys.g = [2.008 2.0061 2.0027];
Sys.Nucs = '14N';
Sys.A = [16 16 86];
Sys.lw = 0.1;
Exp.mwFreq = 9.8;

tc = 10.^(-9:0.5:-7);

for q = 1:numel(tc)
  Sys.tcorr = tc(q);
  [x,y(q,:)] = chili(Sys,Exp);
  y(q,:) = y(q,:)/max(y(q,:));
end

data.x = x;
data.y = y;

if ~isempty(olddata)
  if opt.Display
    subplot(3,1,[1 2]);
    plot(data.x,data.y,'r',olddata.x,olddata.y,'g');
    axis tight;
    subplot(3,1,3);
    plot(data.x,data.y-olddata.y);
  end
  ok = areequal(y,olddata.y,1e-4,'abs');
  err = ~ok;
else
  err = [];
end
