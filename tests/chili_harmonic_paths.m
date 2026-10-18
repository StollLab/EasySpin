function ok = test(opt)

% Check that the first and second harmonic equal the derivatives of the
% absorption spectrum in chili, with and without Gaussian convolution.

Sys.g = [2.008 2.006 2.003];
Sys.tcorr = 1e-9;
Exp.mwFreq = 9.5;
Exp.Range = [334 342];
Exp.nPoints = 1000;

lws = {[0 0], [0.3 0]};

for k = numel(lws):-1:1
  Sys.lw = lws{k};
  for h = 0:2
    Exp.Harmonic = h;
    [x,spc{h+1}] = chili(Sys,Exp);
  end
  d1 = numderiv(x,spc{1});
  d2 = numderiv(x,d1);
  idx = round(0.05*numel(x)):round(0.95*numel(x));
  ok(k) = areequal(spc{2}(idx),d1(idx),1e-2,'rel') && ...
          areequal(spc{3}(idx),d2(idx),1e-2,'rel');

  if opt.Display
    subplot(numel(lws),2,2*k-1);
    plot(x,spc{2},x,d1,'--');
    subplot(numel(lws),2,2*k);
    plot(x,spc{3},x,d2,'--');
  end
end

end

function dy = numderiv(x,y)
d = diff(y)/(x(2)-x(1));
dy = (d([1 1:end]) + d([1:end end]))/2;
end
