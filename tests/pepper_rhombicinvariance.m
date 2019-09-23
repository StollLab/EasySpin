function [err,data] = test(opt,olddata)

% Rhombic symmetry invariance

Sys = struct('S',.5,'g',[1.9,2,2.3]);
Exp = struct('Range',[285 365],'mwFreq',9.5,'nPoints',1054,'Harmonic',0);
Opt = struct('nKnots',[20 2]);
Symmetry = {'Ci','D2h','C2h'};
Sys.lw = 1;

for k = 1:length(Symmetry)
  Opt.Symmetry = Symmetry{k};
  [x,y(k,:)] = pepper(Sys,Exp,Opt);    
end

y = y/max(y(:));
thr = 0.005;
err = ~areequal(y(1,:),y(2,:),thr,'abs') | ~areequal(y(2,:),y(3,:),thr,'abs');

if (opt.Display)
  plot(x,y);
  title('Test 9: Rhombic symmetry invariance');
  legend(Symmetry{:});
  xlabel('magnetic field [mT]');
end

data = [];
