function [err,data] = saffron_hyscore_ganiso(opt,olddata)

Sys.g = [2 6];
Sys.Nucs = '1H';
Sys.A = [2 25];

Exp.Sequence = 'HYSCORE';
Exp.Field = 330;
Exp.tau = 0.001;
Exp.dt = 0.01;
Exp.nPoints = 200;

Opt.nKnots = 181;
Opt.Verbosity = 0;

y = saffron(Sys,Exp,Opt);
y = y/max(abs(y(:)));
y = real(y);

if (opt.Display)
  xlabel('time [us]');
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
  err = [];
else
  ok = areequal(olddata.y,y,1e-10,'rel');
  err = ~ok;
end
