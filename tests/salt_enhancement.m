function [ok,data] = test(opt,olddata)

% Options: Intensity, Hyperfine enhancement
  
Sy = struct('Nucs','14N','g',[2.25 2.25 2],...
  'A',[44,49,54],'Q',[1 1 -2]*0.84,'HStrain',[1 1 1]*100,'lwEndor',0);
Ex = struct('Range',[16 65],'Field',326.5,'nPoints',2001,'mwFreq',9.7);
Si = struct('Threshold',1e-4,'nKnots',[31 12],'Verbosity',opt.Verbosity);

Si.Intensity = 'on';
Si.Enhancement = 'on';
[a,b1] = salt(Sy,Ex,Si);

Si.Intensity = 'on';
Si.Enhancement = 'off';
[a,b2] = salt(Sy,Ex,Si);

Si.Intensity = 'off';
Si.Enhancement = 'off';
[a,b3] = salt(Sy,Ex,Si);

if (opt.Display)
  subplot(3,1,[1 2]);
  plot(a,renorm(b1),'b',a,renorm(b2),'r',a,renorm(b3),'g');
  axis tight
  title('Intensity computations');
  legend('Int+Enh','Int only','all off');
  xlabel('frequency [MHz]');
  if ~isempty(olddata)
    subplot(3,1,3);
    plot(a,b1-olddata.b1,'b',a,b2-olddata.b2,'r',a,b3-olddata.b3,'g');
  end
end

data.b1 = b1;
data.b2 = b2;
data.b3 = b3;

if isempty(olddata)
  ok = [];
else
  e = 1e-6;
  ok = ...
    areequal(b1,olddata.b1,e,'rel') && ...
    areequal(b2,olddata.b2,e,'rel') && ...
    areequal(b3,olddata.b3,e,'rel');
end
