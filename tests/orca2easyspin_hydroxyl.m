function ok = test(opt)

% Read binary ORCA property files

BaseDir = 'orca/';

Files{1} = 'hydroxyl_g.oof';
Files{2} = 'hydroxyl_gA.oof';
Files{3} = 'hydroxyl_gAiso.oof';
Files{4} = 'hydroxyl_gQ.oof';
Files{5} = 'hydroxyl_Q.oof';
Files{6} = 'hydroxyl_HO.oof';

hfCutoff = 1;

for iFile = 1:numel(Files)
  thisFile = Files{iFile};
  if (opt.Display)
    disp(thisFile);
  end
  clear data pars;
  thisFileName = [BaseDir thisFile];
 
  Sys = orca2easyspin(thisFileName);
  Sys = orca2easyspin(thisFileName,hfCutoff);
  
  readerr(iFile) = false;
  if ~isstruct(Sys)
    readerr(iFile) = true;
    disp([  '   error reading ' Files(iFile).name]);
  end
  if (opt.Display)
    Sys
    pause
  end
end

ok = ~readerr;
