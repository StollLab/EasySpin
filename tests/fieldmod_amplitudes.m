function [err,data] = test(opt,olddata)

%=======================================================
% Test 1: various modulation amplitudes
%=======================================================
x = linspace(-10,10,1e3);
y = lorentzian(x,0,1);
ModAmpl = [.05 .1 .2 .5 1:10];
for k = 1:numel(ModAmpl)
  yy(k,:) = fieldmod(x,y,ModAmpl(k),1);
end

if (opt.Display)
  plot(x,yy);
  pause;
end

data.yy = yy;

if ~isempty(olddata)
  ok = areequal(yy,olddata.yy);
  err = ~ok;
else
  err = [];
end
