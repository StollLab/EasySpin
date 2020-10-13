function [ok,data] = test(opt,olddata)

% Regression test for Gaussian lineshape
%-----------------------------------------------------

x = linspace(-1,2,1001);
x0 = 0.9;
fwhm = 0.2;
[ya0,yd0] = gaussian(x,x0,fwhm,0);
[ya1,yd1] = gaussian(x,x0,fwhm,0);

data.ya0 = ya0;
data.yd0 = yd0;
data.ya1 = ya1;
data.yd1 = yd1;

% Plotting
if opt.Display
  subplot(2,1,1)
  plot(x,ya0,x,yd0);
  subplot(2,1,2)
  plot(x,ya1,x,yd1);
end

% Comparison
if ~isempty(olddata)
  thr = 1e-10;
  ok(1) = areequal(data.ya0,olddata.ya0,thr,'rel');
  ok(2) = areequal(data.yd0,olddata.yd0,thr,'rel');
  ok(3) = areequal(data.ya1,olddata.ya1,thr,'rel');
  ok(4) = areequal(data.yd1,olddata.yd1,thr,'rel');
else
  ok = [];
end
