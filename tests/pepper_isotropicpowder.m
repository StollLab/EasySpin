function ok = test(opt)

% Simulation of isotropic powder spectrum

Sys = struct('S',1/2,'g',[2 2 2],'Nucs','63Cu','A',[40 40 40],'HStrain',[1 1 1]*10);
Exp = struct('mwFreq',9.7979,'Range',[340 360],'nPoints',10000);
Opt = struct('Verbosity',0,'GridSize',[30 1]);

Opt.Output = 'summed';
[x,y1] = pepper(Sys,Exp,Opt);
Opt.Output = 'separate';
[x,y2] = pepper(Sys,Exp,Opt);

ok = (size(y1,1)==1) && (size(y2,1)==4) && areequal(y1/max(y1),sum(y2)/max(y1),1e-4,'abs');

if opt.Display
  subplot(3,1,1);
  plot(x,y1);
  xlabel('magnetic field [mT]');
  title('Isotropic powder spectrum, summed output');
  subplot(3,1,2);
  plot(x,y2);
  xlabel('magnetic field [mT]');
  title('Isotropic powder spectrum, separate output');
  subplot(3,1,3);
  plot(x,y1-sum(y2));
  xlabel('magnetic field [mT]');
  title('Difference');
end

