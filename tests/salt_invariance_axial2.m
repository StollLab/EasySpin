function [ok,data] = test(opt,olddata)

% Axial system, symmetry invariance

Sy.S = 1/2;
Sy.Nucs = '1H';
Sy.g = 2;
Sy.A = 2*[-1 -1 2];
Sy.lwEndor = 0.05;

Ex.Range = [10 20];
Ex.Field = 326.5;

Op.GridSize = [30 3];
Op.Verbosity = opt.Verbosity;

%Sy = struct('S',1/2,'Nucs','1H','g',[2 2 2],'A',2*[-1,-1,2],'lwEndor',0.05);
%Ex = struct('Range',[10 20],'Field',326.5,'mwFreq',9.7);
%Op = struct('Threshold',1e-5,'GridSize',[30 5],'Intensity','on',...
%  'Enhancement','off','Verbosity',opt.Verbosity);
Symmetry = {'Dinfh','D6h','D4h','Oh','D3d','Th','D2h',...
    'C4h','C6h','S6','C2h','Ci'};

for k = 1:length(Symmetry)
  Op.GridSymmetry = Symmetry{k};
  [x,y(k,:)] = salt(Sy,Ex,Op);
end

if (opt.Display)
  cla
  plot(x,y);
  title('Axial system, symmetry invariance');
  legend(Symmetry{:});
  
  if ~isempty(olddata)
    hold on
    plot(x,olddata.y,'.k');
  end
end

data.y = y;

if isempty(olddata)
  ok = [];
else
  % Direct test
  ymean = mean(y);
  s = max(ymean);
  ymean = ymean/s;
  yy = y/s;
  dy = (yy - repmat(ymean,size(y,1),1)).^2;
  rms = sqrt(mean(dy,2));
  ok = all(rms<0.01);
  % Regression test
  ok = ok & areequal(olddata.y,y,1e-5,'abs');
end
