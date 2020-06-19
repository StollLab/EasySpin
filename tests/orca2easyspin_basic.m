function ok = test(opt)

% Read binary ORCA property files

BaseDir = 'orca/';

Files{1} = 'tp051512_opt';
Files{2} = 'dioxygen_g';
Files{3} = 'nx4me';

hfCutoff = 1; % MHz

for iFile = 1:numel(Files)
  thisFile = Files{iFile};
  if (opt.Display)
    disp(thisFile);
  end
  clear data pars;
  thisFileName = [BaseDir thisFile];
  
  [~] = orca2easyspin(thisFileName);
  Sys = orca2easyspin(thisFileName,hfCutoff);
  
  readerr(iFile) = false;
  if ~isstruct(Sys)
    readerr(iFile) = true;
    disp([  '   error reading ' Files(iFile).name]);
  end
  if opt.Display
    Sys
    pause
  end
end

ok = ~readerr;
