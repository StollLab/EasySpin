function [err,data] = test(opt,olddata)

% Temperature dependence

Sys = struct('S',1,'g',[1 1 1]*2,'D',200*[1 1 -2],'lw',1);
Exp = struct('Range',[300 380],'mwFreq',9.5,'Harmonic',0);
Opt = struct('nKnots',20);

Temps = [20 10 5 2 1 0.5 0.2];

for t = 1:numel(Temps)
  Exp.Temperature = Temps(t);
  [x,y(t,:)] = pepper(Sys,Exp,Opt);
  leg{t} = sprintf('%g K',Temps(t));
end

if (opt.Display)
  plot(x,y);
  legend(leg{:});
  title('Temperature dependence');
  xlabel('magnetic field [mT]');
  if ~isempty(olddata)
    figure
    plot(x,y-olddata.y);
    title('Residuals');
  end
end

data.y = y;

if ~isempty(olddata)
  err = ~areequal(y,olddata.y,1e-3,'rel');
else
  err = [];
end
