function ok = test(opt)

% Rhombic symmetry invariance

Sys = struct('S',0.5,'g',[1.9,2,2.3]);
Exp = struct('Range',[285 365],'mwFreq',9.5,'nPoints',1054,'Harmonic',0);
Opt = struct('GridSize',[20 0]);
Symmetry = {'D2h','Ci','C1','C2h'};
Sys.lw = 1;

for k = 1:length(Symmetry)
  Opt.GridSymmetry = Symmetry{k};
  [B,spc(k,:)] = pepper(Sys,Exp,Opt);    
end

spc = spc/max(spc(:));
thr = 0.06;
ok(1) = areequal(spc(2,:),spc(1,:),thr,'abs');
ok(2) = areequal(spc(3,:),spc(1,:),thr,'abs');
ok(3) = areequal(spc(4,:),spc(1,:),thr,'abs');

if opt.Display
  subplot(3,1,[1 2]);
  plot(B,spc);
  title('Rhombic symmetry invariance');
  legend(Symmetry{:});
  xlabel('magnetic field [mT]');
  subplot(3,1,3);
  plot(B,spc(2,:)-spc(1,:),B,spc(3,:)-spc(1,:),B,spc(4,:)-spc(1,:));
  legend(Symmetry{2:end});
end
