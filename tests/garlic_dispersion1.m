function [err,data] = test(opt,olddata)

% garlic: dispersion spectrum, first harmonic

Sys = struct('g',2,'Nucs','1H','lw',[0,0.03],'A',10);
Exp = struct('mwFreq',9.7);
Exp.Range = [346 347];
Exp.Harmonic = 1;

Exp.mwPhase = pi/2;
[x,y] = garlic(Sys,Exp);

data.x = x;
data.y = y;

if ~isempty(olddata)
  err = ~areequal(x,olddata.x,1e-10,'abs') || ~areequal(y,olddata.y,1e-4,'rel');
else
  err = [];
end

if (opt.Display)
  if ~isempty(olddata)
    plot(olddata.x,olddata.y,'b',data.x,data.y,'r');
    legend('old','new');
  else
    plot(x,y);
  end
  xlabel('magnetic field (mT)');
  axis tight;
end
