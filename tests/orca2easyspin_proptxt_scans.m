function ok = test()

% Read text-based property file for ORCA 5.0.x

BaseDir = 'orca/scans';

files{1} = 'hydroxyl_optscan_property.txt';
files{2} = 'hydroxyl_parascan_property.txt';

for k = 1:numel(files)
  filename = [BaseDir filesep files{k}];
  Sys = orca2easyspin(filename);
  ok(k) = true;
end
