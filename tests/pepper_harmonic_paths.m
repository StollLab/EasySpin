function ok = test(opt)

% Check that the first and second harmonic equal the derivatives of the
% absorption spectrum, independent of how pepper constructs the spectrum
% (template summation, projection, anisotropic widths, field and frequency
% sweeps).

Exp.mwFreq = 9.5;
Exp.Range = [315 345];
Exp.nPoints = 2000;

SysIso.g = 2.02;
SysPowder.g = [2 2.05 2.1];

% {Sys, lw, lwAnisoG, Exp}
cases = {
  SysIso,    [0.5 0],   [],           Exp  % summation, Gaussian template
  SysIso,    [0 0.5],   [],           Exp  % summation, Lorentzian template
  SysIso,    [0.5 0.3], [],           Exp  % summation, Gaussian template + Lorentzian convolution
  SysPowder, [1 0],     [],           Exp  % projection, Gaussian convolution
  SysPowder, [0 1],     [],           Exp  % projection, Lorentzian convolution
  SysPowder, [1 1],     [],           Exp  % projection, both convolutions
  SysPowder, [0 0],     [50 100 150], Exp  % anisotropic widths only
  };

% Frequency sweep
ExpF.Field = 340;
ExpF.mwRange = [9.3 9.8];
ExpF.nPoints = 2000;
cases(end+1,:) = {SysPowder, [0 20], [], ExpF};

for k = size(cases,1):-1:1
  Sys = cases{k,1};
  Sys.lw = cases{k,2};
  if ~isempty(cases{k,3}), Sys.lwAnisoG = cases{k,3}; end
  Exp_ = cases{k,4};
  for h = 0:2
    Exp_.Harmonic = h;
    [x,spc{h+1}] = pepper(Sys,Exp_);
  end
  d1 = numderiv(x,spc{1});
  d2 = numderiv(x,d1);
  % Exclude the range edges, where convolution and numerical
  % differentiation have different artifacts
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
