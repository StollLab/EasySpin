function ok = test()

% Legacy ORCA 3.0.3 files: check that main output files (.oof) and binary
% property files (.prop) from the same calculation give the same tensors.

folder = 'orca/v3.0.3/';

calcs = {'dioxygen_g','dioxygen_gD','hydroxyl_g','hydroxyl_gA', ...
  'hydroxyl_gAiso','hydroxyl_gQ','hydroxyl_Q','hydroxyl_HO','hydroxyl_098_v303'};

T = @(v,ang) erot(ang).'*diag(v)*erot(ang);
maxdiff = @(X,Y) max(abs(X(:)-Y(:)));

for k = numel(calcs):-1:1
  SysM = orca2easyspin([folder calcs{k} '.oof']);
  SysP = orca2easyspin([folder calcs{k} '.prop']);

  okk = true;
  for f = {'g','D','A','Q'}
    okk = okk && isfield(SysM,f{1})==isfield(SysP,f{1});
  end
  if isfield(SysM,'g')
    okk = okk && maxdiff(T(SysM.g,SysM.gFrame),T(SysP.g,SysP.gFrame))<1e-6;
  end
  if isfield(SysM,'D')
    okk = okk && maxdiff(T(SysM.D,SysM.DFrame),T(SysP.D,SysP.DFrame))<1;  % MHz
  end
  if isfield(SysM,'A')
    okk = okk && maxdiff(T(SysM.A,SysM.AFrame),T(SysP.A,SysP.AFrame))<0.01;  % MHz
  end
  if isfield(SysM,'Q')
    okk = okk && maxdiff(T(SysM.Q,SysM.QFrame),T(SysP.Q,SysP.QFrame))<1e-4;  % MHz
  end

  ok(k) = okk;
end
