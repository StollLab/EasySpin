function [ok,data] = test(opt,olddata)

Sys.Nucs = '14N';
%Sys.A = [-1 -1 2]*2+1.3; Sys.Q = [-1 -1 2]*0.3;
Sys.A = [-1 -1 2]*0.5+0.8;
Sys.Q = [-1 -1 2]*0.2;

Exp.Sequence = 'HYSCORE';
Exp.Field = 330;
Exp.tau = 0.001;
Exp.dt = 0.1;
Exp.nPoints = 200;

Opt.GridSize = 181;
Opt.Expand = 3;

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
  end
end

data.y = y;

if isempty(olddata)
  ok = [];
else
  ok = areequal(olddata.y,y,1e-10,'rel');
end
