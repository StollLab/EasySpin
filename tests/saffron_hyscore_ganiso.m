function [ok,data] = saffron_hyscore_ganiso(opt,olddata)

% Regression test: powder HYSCORE of one 1H coupled to an electron spin with
% strongly anisotropic g. The g anisotropy tilts the electron spin
% quantization axis away from the magnetic field, which changes the
% hyperfine field at the nucleus. Without Exp.mwFreq, all orientations are
% excited. The real part of the normalized signal is compared, so the test
% checks the shape of the signal, not its absolute amplitude.

Sys.g = [2 6];  % axial, unrealistically anisotropic
Sys.Nucs = '1H';
Sys.A = [2 25];  % MHz

Exp.Sequence = 'HYSCORE';
Exp.Field = 330;  % mT
Exp.tau = 0.001;  % µs, short enough to avoid blind spots
Exp.dt = 0.01;  % µs
Exp.nPoints = 200;

Opt.GridSize = 181;
Opt.Verbosity = 0;

y = saffron(Sys,Exp,Opt);
y = y/max(abs(y(:)));
y = real(y);

if (opt.Display)
  if ~isempty(olddata)
    figure;
    subplot(1,2,1);
    pcolor(y); shading flat; title('new');
    subplot(1,2,2);
    pcolor(olddata.y); shading flat; title('old');
  else
    figure;
    pcolor(y); shading flat; title('new');    
  end
end

data.y = y;

if isempty(olddata)
  ok = [];
else
  ok = areequal(olddata.y,y,1e-10,'rel');
end
