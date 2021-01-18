function [ok,data] = test(opt,olddata)

% Rotational invariance
      
Sy = struct('S',1/2,'Nucs','1H','g',[2 2.01 2.02],...
  'A',2*[-1,-1,2],'lwEndor',0.1);
Ex = struct('Range',[10 20],'Field',326.5,'mwFreq',9.7);
Si = struct('Threshold',1e-5,'GridSize',[20 5],'Intensity','on',...
  'Enhancement','off','Verbosity',opt.Verbosity);

[a,y(1,:)] = salt(Sy,Ex,Si);
Sy.AFrame = pi/180*[23 72 -12];
[a,y(2,:)] = salt(Sy,Ex,Si);
Sy.AFrame = pi/180*[0 9 65];
[a,y(3,:)] = salt(Sy,Ex,Si);

if (opt.Display)
  plot(a,y);
  title('Rotation invariance, axial A');
  legend('axial A along xyz (D2h)','A tilted (Ci)','A tilted again (Ci)');
  xlabel('frequency [MHz]');
end

data.y = y;

if isempty(olddata)
  ok = [];
else
  % Direct test for invariance
  ymean = mean(y);
  s = max(ymean);
  ymean = ymean/s;
  yy = y/s;
  dy = (yy - repmat(ymean,size(yy,1),1)).^2;
  rms = sqrt(mean(dy,2));
  ok = all(rms<0.01);
  % Regression test
  ok = ok && areequal(data.y,y,1e-10,'rel');
end
