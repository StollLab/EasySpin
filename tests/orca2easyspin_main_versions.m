function ok = test(opt)

% Read ORCA output files

BaseDir = 'orca/';

Files{1} = 'hydroxyl_098_v303.oof';
Files{2} = 'hydroxyl_098_v400.oof';
Files{3} = 'hydroxyl_098_v503.oof';

for k = 1:numel(Files)
  filename = [BaseDir filesep Files{k}];
  Sys = orca2easyspin(filename);
  ok(k) = true;
end
