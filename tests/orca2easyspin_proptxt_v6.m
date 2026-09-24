function ok = test()

% Compare spin systems read from ORCA 6 property files (*.property.txt)
% with those read from the corresponding main output files.

names = {'nitroxide_tilted_v610','tripletformaldehyde_tilted_v610'};
folder = 'orca/';

tensor = @(v,ang) erot(ang).'*diag(v)*erot(ang);
maxdiff = @(T1,T2) max(abs(T1(:)-T2(:)));

for k = numel(names):-1:1
  SysP = orca2easyspin([folder names{k} '.property.txt']);
  SysM = orca2easyspin([folder names{k} '.out']);

  okk = SysP.S==SysM.S && SysP.charge==0;
  okk = okk && maxdiff(SysP.xyz,SysM.xyz)<1e-5;
  okk = okk && maxdiff(tensor(SysP.g,SysP.gFrame),tensor(SysM.g,SysM.gFrame))<1e-6;
  okk = okk && isfield(SysP,'D')==isfield(SysM,'D');
  if isfield(SysP,'D')
    okk = okk && maxdiff(tensor(SysP.D,SysP.DFrame),tensor(SysM.D,SysM.DFrame))<1;  % MHz
  end
  okk = okk && strcmp(SysP.Nucs,SysM.Nucs) && isequal(SysP.NucsIdx,SysM.NucsIdx);
  for n = 1:numel(SysP.NucsIdx)
    okk = okk && maxdiff(tensor(SysP.A(n,:),SysP.AFrame(n,:)),tensor(SysM.A(n,:),SysM.AFrame(n,:)))<0.05;  % MHz
    okk = okk && maxdiff(tensor(SysP.Q(n,:),SysP.QFrame(n,:)),tensor(SysM.Q(n,:),SysM.QFrame(n,:)))<1e-4;  % MHz
  end
  ok(k) = okk;
end

end
