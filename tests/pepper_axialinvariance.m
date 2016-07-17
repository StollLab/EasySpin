function [err,data] = test(opt,olddata)

% An axial spectrum should be independent of
% the region of orientational integration.

Sys = struct('S',.5,'g',[1.9,2.3]);
Exp = struct('Range',[285 365],'mwFreq',9.5,'nPoints',1024,'Harmonic',0);
Opt = struct('nKnots',[19 5],'Method','perturb1');
Symmetry = {'Dinfh','D6h','D4h','Oh','D3d','Th','D2h',...
    'C4h','C6h','S6','C2h','Ci'};
Sys.lw = 1;

for k = 1:numel(Symmetry)
  Opt.Symmetry = Symmetry{k};
  [x,y(k,:)] = pepper(Sys,Exp,Opt);
  dy(k) = max(y(k,:) - y(1,:));
end

y = y/max(y(:));

if (opt.Display)
  plot(x,y);
  title('Axial symmetry invariance');
  xlabel('magnetic field [mT]');
  ylabel('normalized amplitude');
end

for k = 1:numel(Symmetry)
  dy(k) = max(abs(y(k,:) - y(end,:)));
end

err = max(dy)>0.005;

data = [];
