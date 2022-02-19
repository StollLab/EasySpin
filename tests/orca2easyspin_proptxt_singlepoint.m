function ok = test()

% Read text-based property file for ORCA 5.0.x

BaseDir = 'orca';

files{1} = 'water_property.txt';
files{2} = 'dioxygen_singlepoint_property.txt';
files{3} = 'methyl_all_property.txt';

S = [0 1 1/2];

for k = 1:numel(files)
  filename = [BaseDir filesep files{k}]
  Sys = orca2easyspin(filename);
  ok(k) = Sys.S==S(k);
end
