function [ok,data] = test(opt,olddata)

%=======================================================
% Test various modulation harmonics
%=======================================================
x = linspace(-10,10,1e3);
y = gaussian(x,-2,1);
for Harmonic = 1:4
  yy(Harmonic,:) = fieldmod(x,y,1,Harmonic);
end

if (opt.Display)
  plot(x,yy);
  pause;
end

data.yy = yy;

if ~isempty(olddata)
  ok = areequal(yy,olddata.yy,1e-10,'rel');
else
  ok = [];
end
