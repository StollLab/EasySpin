function ok = test()

% Read text-based property file for ORCA 5.x

BaseDir = 'orca/5.0.2';

Sys = orca2easyspin([BaseDir filesep 'dioxygen_g.oof']);
Sys = orca2easyspin([BaseDir filesep 'hydroxyl_gA.oof']);

ok = true;
