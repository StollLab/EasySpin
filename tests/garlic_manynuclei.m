function [err,data] = test(opt,olddata)

%=======================================================
% Many nuclei
%=======================================================
System = struct('g',2,'Nucs','1H,14N','A',[30,40],'n',[5 4]);
System.lw = [0 0.1]; % only Lorentzian broadening
Exper = struct('mwFreq',9.7,'nPoints',2000);
[x,y] = garlic(System,Exper);

data.x = x;
data.y = y;

if ~isempty(olddata)
  ok = areequal(x,olddata.x,1e-10,'rel') && areequal(y,olddata.y,1e-10,'rel');
  err = ~ok;
else
  err = [];
end

if (opt.Display)
  if ~isempty(olddata)
    subplot(3,1,[1 2]);
    plot(olddata.x,olddata.y,'g',data.x,data.y,'r');
    legend('old','new');
    subplot(3,1,3);
    plot(data.x,olddata.y-data.y);
  else
    plot(x,y);
  end
  xlabel('magnetic field [mT]');
  axis tight;
end
