function ok = test(opt)

% Check that the first and second harmonic equal the derivatives of the
% absorption spectrum, independent of the accumulation method in garlic
% and in the fast-motion regime.

Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = 20;
Exp.mwFreq = 9.5;
Exp.Range = [334 345];
Exp.nPoints = 2000;

% {lw, tcorr, AccumMethod}
% Line widths are large compared to the increment, so that the numerical
% derivatives are accurate.
cases = {
  [0.2 0.15], [],    'explicit'
  [0 0.3],    [],    'template'
  [0.2 0.15], [],    'binning'
  [0.2 0],    [],    'binning'
  [0 0.3],    1e-10, 'explicit'  % fast motion
  };
Sys.A = [10 20 30];  % only used in the fast-motion case

for k = size(cases,1):-1:1
  Sys_ = Sys;
  Sys_.lw = cases{k,1};
  if isempty(cases{k,2})
    Sys_.A = 20;
  else
    Sys_.tcorr = cases{k,2};
  end
  Opt.AccumMethod = cases{k,3};
  for h = 0:2
    Exp.Harmonic = h;
    [x,spc{h+1}] = garlic(Sys_,Exp,Opt);
  end
  d1 = numderiv(x,spc{1});
  d2 = numderiv(x,d1);
  idx = round(0.05*numel(x)):round(0.95*numel(x));
  ok(k) = areequal(spc{2}(idx),d1(idx),1e-2,'rel') && ...
          areequal(spc{3}(idx),d2(idx),1e-2,'rel');

  if opt.Display
    subplot(size(cases,1),2,2*k-1);
    plot(x,spc{2},x,d1,'--');
    subplot(size(cases,1),2,2*k);
    plot(x,spc{3},x,d2,'--');
  end
end

end

function dy = numderiv(x,y)
d = diff(y)/(x(2)-x(1));
dy = (d([1 1:end]) + d([1:end end]))/2;
end
