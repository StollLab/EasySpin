function ok = test(opt)

% Make sure the 'lsq' mode works correctly.

rng(232)
x = linspace(0,2,1001);
y0 = gaussian(x,0.7,0.5);
y = addnoise(y0,5,'n');
scale0 = 2;
y = y/scale0;

yref = y0;

[ys,scale_fit] = rescaledata(y,yref,'lsq');

scale = linspace(0.8,1.2,10001)*scale0;
for i = numel(scale):-1:1
  ssd(i) = sum((y-yref/scale(i)).^2);
end

[~,idx] = min(ssd);
scale_min = scale(idx);

ok = areequal(scale_fit,scale_min,1e-4,'rel');

if opt.Display
  subplot(2,1,1)
  plot(x,ys,x,yref);
  axis tight
  legend('rescaled data','reference');
  subplot(2,1,2)
  plot(scale,ssd);
  xlabel('scale factor');
  ylabel('ssd(yscaled,yref)');
  xline(scale0,'Label','true scale','LineStyle',':');
  xline(scale_fit,'Label','fitted scale','Color','r');
  xline(scale_min,'Label','ssd minimum','LineStyle','--');
end
