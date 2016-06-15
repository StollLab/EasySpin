function [err,data] = test(opt,olddata)

% Excitation profile checks
  
Sy = struct('S',1/2,'Nucs','14N','g',[2.2 2.2 2],...
  'A',[4,4,5],'Q',-0.84*[-1,-1,2],'HStrain',[1 1 1]*1,'lwEndor',0.1);
Ex = struct('Range',[0 35],'Field',3394,'nPoints',2^12,'mwFreq',95);
Si = struct('Threshold',1e-4,'nKnots',50,'Intensity','on');
Si.Verbosity = opt.Verbosity;

[a,b1] = salt(Sy,Ex,Si);
Ex.ExciteWidth = 3e3;
[a,b2] = salt(Sy,Ex,Si);
Ex.mwFreq = 104.51;
[a,b3] = salt(Sy,Ex,Si);

if (opt.Display)
  subplot(3,1,[1 2]);
  h = plot(a,renorm(b1),'b',a,renorm(b2),'r',a,renorm(b3),'g');
  set(h(1),'LineWidth',2);
  legend('full excitation','300MHz at 95Ghz (gz)','300MHz at 104.51GHz (gx)');
  xlabel('frequency [MHz]');
  if ~isempty(olddata)
    subplot(3,1,3);
    plot(a,renorm(b1)-renorm(olddata.b1),'b',a,renorm(b2)-renorm(olddata.b2),'r',...
      a,renorm(b3)-renorm(olddata.b3),'g');
  end
end

data.b1 = b1;
data.b2 = b2;
data.b3 = b3;

if isempty(olddata)
  err = [];
else
  e = max([b1 b2 b3])*1e-6;
  ok = ...
    areequal(b1,olddata.b1,e) & ...
    areequal(b2,olddata.b2,e) & ...
    areequal(b3,olddata.b3,e);
  err = ~ok;
end
