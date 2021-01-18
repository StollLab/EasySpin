function ok = test(opt)

% Rhombic symmetry invariance

Sys = struct('S',.5,'g',[1.9,2,2.3]);
Exp = struct('Range',[285 365],'mwFreq',9.5,'nPoints',1054,'Harmonic',0);
Opt = struct('GridSize',[20 2]);
Symmetry = {'D2h','Ci','C1','C2h'};
Sys.lw = 1;

for k = 1:length(Symmetry)
  Opt.GridSymmetry = Symmetry{k};
  [x,y(k,:)] = pepper(Sys,Exp,Opt);    
end

y = y/max(y(:));
thr = 0.005;
ok(1) = areequal(y(2,:),y(1,:),thr,'abs');
ok(2) = areequal(y(3,:),y(1,:),thr,'abs');
ok(3) = areequal(y(4,:),y(1,:),thr,'abs');

if opt.Display
  plot(x,y);
  title('Test 9: Rhombic symmetry invariance');
  legend(Symmetry{:});
  xlabel('magnetic field [mT]');
end
